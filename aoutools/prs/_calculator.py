""""PRS calculator"""

from typing import Optional
import logging
import hail as hl
import hailtop.fs as hfs
import pandas as pd
from aoutools._utils.helpers import SimpleTimer
from ._utils import _log_timing
from ._calculator_utils import (
    _prepare_samples_to_keep,
    _orient_weights_for_split,
    _check_allele_match,
    _calculate_dosage,
    _prepare_weights_for_chunking,
    _create_1bp_intervals,
)
from ._config import PRSConfig

logger = logging.getLogger(__name__)


def _prepare_mt_split(
    vds: hl.vds.VariantDataset,
    weights_table: hl.Table,
    config: PRSConfig
) -> hl.MatrixTable:
    """
    Prepares a MatrixTable for the split-multi PRS calculation path.


    Prepares a MatrixTable for split-multi PRS calculation.

    Splits multi-allelic sites in the VDS, orients weights based on
    `ref_is_effect_allele`, joins with the weights table using (locus, alleles),
    and calculates dosage using GT after splitting.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        An interval-filtered VariantDataset.
    weights_table : hail.Table
        A chunk of the PRS weights table.
    config : PRSConfig
        A configuration object controlling PRS behavior, including
        `split_multi`, `ref_is_effect_allele`, and `detailed_timings`.

    Returns
    -------
    hail.MatrixTable
        A MatrixTable annotated with `weights_info` and per-variant `dosage`.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    with _log_timing(
        "Planning: Splitting multi-allelic variants and joining",
        config.detailed_timings
    ):
        mt = hl.vds.split_multi(vds).variant_data

        weights_ht_processed = _orient_weights_for_split(weights_table, config)
        mt = mt.annotate_rows(
            weights_info=weights_ht_processed[mt.row_key]
        )

        mt = mt.filter_rows(hl.is_defined(mt.weights_info))

    with _log_timing(
        "Planning: Calculating per-variant dosage",
        config.detailed_timings,
    ):
        # After splitting, LGT is converted to GT, so we can
        # directly and safely use the built-in dosage calculator.
        # See the source code for `hl.vds.split_multi` for details.
        mt = mt.annotate_entries(dosage=mt.GT.n_alt_alleles())

        return mt


def _prepare_mt_non_split(
    vds: hl.vds.VariantDataset,
    weights_table: hl.Table,
    config: PRSConfig
) -> hl.MatrixTable:
    """
    Prepares a MatrixTable for the non-split PRS calculation path.

    This function takes an interval-filtered VDS and joins it with the
    weights table using a locus-based key. It optionally performs a strict
    allele match to handle allele orientation and then calculates dosage
    using a custom multi-allelic dosage function.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        An interval-filtered Variant Dataset.
    weights_table : hail.Table
        A chunk of the weights table.
    config : PRSConfig
        A configuration object controlling `strict_allele_match` and
        `detailed_timings`.

    Returns
    -------
    hail.MatrixTable
        A MatrixTable, filtered and annotated with `weights_info` and `dosage`
        for the specified effect allele.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    mt = vds.variant_data

    with _log_timing(
        "Planning: Annotating variants with weights", config.detailed_timings
    ):
        mt = mt.annotate_rows(weights_info=weights_table[mt.locus])
        mt = mt.filter_rows(hl.is_defined(mt.weights_info))

    if config.strict_allele_match:
        with _log_timing(
                "Planning: Performing strict allele match",
                config.detailed_timings
        ):
            is_valid_pair = _check_allele_match(mt, mt.weights_info)
            mt = mt.filter_rows(is_valid_pair)

    with _log_timing(
        "Planning: Calculating per-variant dosage",
        config.detailed_timings,
    ):
        mt = mt.annotate_entries(dosage=_calculate_dosage(mt))

    return mt


def _calculate_prs_chunk(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    config: PRSConfig
) -> hl.Table:
    """
    Calculates a Polygenic Risk Score (PRS) for a single chunk of variants.

    This function serves as the core computation step. It prepares the variant
    data depending on whether multi-allelic splitting is enabled, and computes
    the PRS using dosage-weight aggregation.

    Parameters
    ----------
    weights_table : hail.Table
        A pre-filtered chunk of the full weights table, keyed by 'locus'.
    vds : hail.vds.VariantDataset
        The Variant Dataset containing genotypes to score.
    config : PRSConfig
        A configuration object specifying settings such as `split_multi`,
        `include_n_matched`, and `sample_id_col`.

    Returns
    -------
    hail.Table
        A Hail Table with one row per sample and a PRS column. If requested,
        also includes the number of matched variants ('n_matched').

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    if config.split_multi:
        mt = _prepare_mt_split(
            vds=vds,
            weights_table=weights_table,
            config=config,
        )
    else:
        mt = _prepare_mt_non_split(
            vds=vds,
            weights_table=weights_table,
            config=config,
        )

    # Chunks aggregation
    prs_table = mt.select_cols(
        prs=hl.agg.sum(mt.dosage * mt.weights_info.weight)
    ).cols()

    if config.include_n_matched:
        with _log_timing(
                "Computing shared variants count", config.detailed_timings
        ):
            # Using hl.agg.count() within the `select_cols` block won't work
            # since homozygous reference are set to missing while `agg.count`
            # counts the number of rows for which that specific sample has a
            # non-missing genotype calls.
            # This is two-pass approach and thus less performant.
            n_matched = mt.count_rows()
            logger.info("%d variants in common in this chunk.", n_matched)
            prs_table = prs_table.annotate(n_matched=n_matched)

    # Rename sample ID column to user-defined name
    prs_table = prs_table.rename({'s': config.sample_id_col})

    # Drop all global annotations to minimize memory footprint
    return prs_table.select_globals()


def _process_chunks(
    full_weights_table: hl.Table,
    n_chunks: int,
    vds: hl.vds.VariantDataset,
    config: PRSConfig
) -> list[pd.DataFrame]:
    """
    Iteratively processes each chunk of the weights table.

    This helper function orchestrates the main PRS calculation loop. For each
    chunk, it filters the Variant Dataset to the relevant genomic intervals,
    computes the PRS using `_calculate_prs_chunk`, and converts the result to a
    Pandas DataFrame.

    Parameters
    ----------
    full_weights_table : hail.Table
        The full weights table, annotated with a 'chunk_id' field.
    n_chunks : int
        The total number of chunks to process.
    vds : hail.vds.VariantDataset
        The Variant Dataset containing genotype data, optionally
        filtered for samples.
    config : PRSConfig
        A configuration object specifying PRS settings, including
        `detailed_timings` and `sample_id_col`.

    Returns
    -------
    list[pd.DataFrame]
        A list of Pandas DataFrames, where each DataFrame contains
        the partial PRS results for one chunk.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    partial_dfs = []
    for i in range(n_chunks):
        with _log_timing(
            f"Processing chunk {i + 1}/{n_chunks}", config.detailed_timings
        ):
            # Use .persist() to avoid recomputation of the same chunk in
            # _calculate_prs_chunk, specifically during:
            # 1. Creation of interval_ht
            # 2. Annotating rows with PRS weight information
            weights_chunk = full_weights_table.filter(
                full_weights_table.chunk_id == i
            ).persist()

        with _log_timing(
            "Planning: Filtering VDS to variants in weights table chunk",
            config.detailed_timings,
        ):
            intervals_to_filter = _create_1bp_intervals(weights_chunk)
            # If filter_intervals filters the main vds and reassigns to vds
            # again, subsequent operation will try to filter empty variable.
            vds_chunk = hl.vds.filter_intervals(
                vds, intervals_to_filter, keep=True
            )

            chunk_prs_table = _calculate_prs_chunk(
                weights_table=weights_chunk,
                vds=vds_chunk,
                config=config
            )

            # Convert the per-chunk Hail Table to a Pandas DataFrame.
            partial_dfs.append(chunk_prs_table.to_pandas())

    return partial_dfs


def _aggregate_and_export(
    partial_dfs: list[pd.DataFrame],
    output_path: str,
    config: PRSConfig
) -> None:
    """
    Aggregates partial Pandas DataFrame results and exports the final result.

    This helper function handles the final aggregation and export stage of the
    PRS pipeline. It concatenates a list of partial DataFrames, groups them by
    sample ID, sums the PRS scores, and writes the final aggregated results to
    a specified cloud storage path.

    Parameters
    ----------
    partial_dfs : list[pd.DataFrame]
        A list of Pandas DataFrames, where each contains partial PRS results
        for a chunk.
    output_path : str
        A destination path on GCS to write the final tab-separated file.
    config : PRSConfig
        A configuration object that specifies `sample_id_col` and
        `detailed_timings`.

    Returns
    -------
    None

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    if not partial_dfs:
        logger.warning(
            "No PRS results were generated. No output file will be created."
        )
        return

    with _log_timing(
            "Aggregating results with Pandas", config.detailed_timings
    ):
        combined_df = pd.concat(partial_dfs, ignore_index=True)
        final_df = combined_df.groupby(config.sample_id_col).sum()

    with _log_timing(
            f"Exporting final result to {output_path}", config.detailed_timings
    ):
        with hfs.open(output_path, 'w') as f:
            final_df.to_csv(f, sep='\t', index=True, header=True)


def calculate_prs(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    output_path: str,
    config: PRSConfig = PRSConfig()
) -> Optional[str]:
    """
    Calculates a Polygenic Risk Score (PRS) and exports the result to a file.

    This function is the main entry point for the PRS calculation workflow. It
    processes a weights table in chunks, using a filter_intervals approach to
    select variants from the VDS for each chunk. Partial results are then
    converted to Pandas DataFrames and aggregated to produce the final score
    file.

    Notes
    -----
    By default (`config.split_multi=True`), this function prioritizes
    robustness over performance by splitting multi-allelic variants.

    This split_multi process includes creating a minimal representation for
    variants. For example, for a variant chr1:10075251 A/G in the weights
    table, split_multi can intelligently match it to a complex indel in the VDS
    (e.g., alleles=['AGGGC', 'A', 'GGGGC']) by simplifying the VDS
    representation to its minimal form (['A', 'G']) for 'AGGGC' -> 'GGGGC'.

    The non-split path (`config.split_multi=False`) is a faster but less robust
    alternative. It relies on a direct string comparison of alleles and will
    fail to match the complex variant described above. Furthermore, if the
    weights table contains multiple entries for the same locus, the non-split
    path will arbitrarily select only one of them. This "power-user" option
    should only be used if you are certain that both your VDS and weights table
    contain only simple, well-matched, bi-allelic variants.

    Parameters
    ----------
    weights_table : hail.Table
        A Hail table containing variant weights. Must contain the following
        columns:
        - 'chr': str
        - 'pos': int32
        - 'effect_allele': str
        - 'noneffect_allele': str
        - A column for the effect weight (float64), specified by
          `weight_col_name`.
    vds : hail.vds.VariantDataset
        A Hail VariantDataset containing both variant and sample data.
    output_path : str
        A GCS path (starting with 'gs://') to write the final tab-separated
        output file.
    config : PRSConfig, optional
        A configuration object for all optional parameters. If not provided,
        default settings will be used. See the `PRSConfig` class for details
        on all available settings.

    Returns
    -------
    str or None
        The output path if results are successfully written; otherwise, None.
        The output file is a tab-separated text file with:
        - A sample ID column (as configured in `config.sample_id_col`)
        - 'prs': The calculated PRS value
        - 'n_matched' (optional): The number of variants used to calculate
          the score, included if `config.include_n_matched` is True.

    Raises
    ------
    ValueError
        If `output_path` is not a valid GCS path, or if the `weights_table`
        is empty after validation.
    TypeError
        If the `config.samples_to_keep` argument is of an unsupported type.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    timer = SimpleTimer()
    with timer:
        if not output_path.startswith('gs://'):
            raise ValueError(
                "The 'output_path' must be a Google Cloud Storage (GCS) "
                "path, starting with 'gs://'."
            )

        logger.info(
            "Starting PRS calculation. Final result will be at: %s",
            output_path,
        )

        if config.samples_to_keep is not None:
            with _log_timing(
                "Planning: Filtering to specified samples",
                config.detailed_timings
            ):
                samples_ht = _prepare_samples_to_keep(config.samples_to_keep)
                vds = hl.vds.filter_samples(vds, samples_ht)

        full_weights_table, n_chunks = _prepare_weights_for_chunking(
            weights_table=weights_table,
            config=config,
            validate_table=True,
        )

        partial_dfs = _process_chunks(
            full_weights_table=full_weights_table,
            n_chunks=n_chunks,
            vds=vds,
            config=config,
        )

        _aggregate_and_export(
            partial_dfs=partial_dfs,
            output_path=output_path,
            config=config,
        )

    # Report the total time using the duration captured by the context manager
    logger.info(
        "PRS calculation complete. Total time: %.2f seconds.", timer.duration
    )
    return output_path if partial_dfs else None
