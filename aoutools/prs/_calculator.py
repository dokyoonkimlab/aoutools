""""PRS calculator"""

import typing
import logging
from math import ceil
import hail as hl
import hailtop.fs as hfs
import pandas as pd
from aoutools._utils.helpers import SimpleTimer
from ._utils import (
    _log_timing,
    _standardize_chromosome_column,
    _prepare_samples_to_keep,
)
from ._config import PRSConfig

logger = logging.getLogger(__name__)


def _validate_and_prepare_weights_table(
    table: hl.Table,
    weight_col_name: str,
    log_transform_weight: bool
) -> hl.Table:
    """
    Validates and prepares a single weights table for PRS calculation.

    This function ensures the table has the required columns with the correct
    types, standardizes the chromosome format, handles the weight column,
    and keys the table by locus for joining with the VDS.

    Parameters
    ----------
    table : hail.Table
        The input weights table. Must contain 'chr', 'pos', 'effect_allele',
        'noneffect_allele', and a weight column.
    weight_col_name : str
        The name of the column containing the effect weights.
    log_transform_weight : bool
        If True, applies a natural log transformation to the weight column.

    Returns
    -------
    hail.Table
        A validated, standardized, and keyed Hail Table.

    Raises
    ------
    TypeError
        If the specified `weight_col_name` does not exist, if other required
        columns are missing, or if any required column has an incorrect
        data type.
    """
    if weight_col_name not in table.row:
        raise TypeError(
            f"Specified weight column '{weight_col_name}' not found in table."
        )
    table = table.rename({weight_col_name: 'weight'})

    required_cols = {
        'chr': hl.tstr,
        'pos': hl.tint32,
        'effect_allele': hl.tstr,
        'noneffect_allele': hl.tstr,
        'weight': hl.tfloat64,
    }
    for col, expected_type in required_cols.items():
        if col not in table.row:
            raise TypeError(f"Weights table is missing required column: '{col}'.")
        if table[col].dtype != expected_type:
            raise TypeError(f"Column '{col}' has incorrect type.")

    table = _standardize_chromosome_column(table)
    if log_transform_weight:
        table = table.annotate(weight=hl.log(table.weight))

    table = table.annotate(
        locus=hl.locus(table.chr, table.pos, reference_genome='GRCh38')
    )
    table = table.key_by('locus')
    return table.select('effect_allele', 'noneffect_allele', 'weight')
    # return table.select('chr', 'pos', 'effect_allele', 'noneffect_allele', 'weight')



def _calculate_dosage(
    mt: hl.MatrixTable,
    score_name: str = ''
) -> hl.expr.Int32Expression:
    """
    Calculates the dosage of the effect allele for a specific score.

    This expression handles both global (GT) and local (LGT/LA) genotype
    encoding formats, which enables sparse storage of homozygous reference
    calls. It correctly computes dosage at multi-allelic sites.

    Parameters
    ----------
    mt : hail.MatrixTable
        MatrixTable annotated with a `weights_info` struct. For batch mode,
        the struct is named `weights_info_{score_name}`.
    score_name : str, optional
        The identifier for the PRS being calculated, used to access the
        correct weights annotation in batch mode.

    Returns
    -------
    hail.expr.Int32Expression
        An expression for the effect allele dosage.
    """
    weights_field_name = (
        f'weights_info_{score_name}' if score_name else 'weights_info'
    )
    effect_allele = mt[weights_field_name].effect_allele
    ref_is_effect = effect_allele == mt.alleles[0]

    # Check for 'GT' field to handle different VDS versions by their
    # genotype encoding scheme.
    if 'GT' in mt.entry:
        # Global-indexed format: 'GT' contains indices that refer
        # directly to the global 'alleles' array.
        # Example: if GT is [0, 1] and mt.alleles is ['A', 'G', 'T'],
        # this expression reconstructs the sample's alleles as ['A', 'G'].
        alleles_expr = hl.or_missing(
            hl.is_defined(mt.GT),
            hl.array([mt.alleles[mt.GT[0]], mt.alleles[mt.GT[1]]])
        )
    else:
        # Local-indexed format: 'LGT' indices refer to the 'LA'
        # (local-to-global) map, which then refers to 'alleles'.
        # Example: LGT=[0, 1], LA=[0, 2], mt.alleles=['A', 'C', 'G']
        # 1. LGT[0] is 0. LA[0] is 0. mt.alleles[0] is 'A'.
        # 2. LGT[1] is 1. LA[1] is 2. mt.alleles[2] is 'G'.
        # The reconstructed alleles are ['A', 'G'].
        alleles_expr = hl.or_missing(
            hl.is_defined(mt.LGT) & hl.is_defined(mt.LA),
            hl.array([
                mt.alleles[hl.or_else(mt.LA[mt.LGT[0]], 0)],
                mt.alleles[hl.or_else(mt.LA[mt.LGT[1]], 0)]
            ])
        )

    # The hl.case statement cleanly handles missing genotypes by assuming
    # they are homozygous reference.
    return hl.case() \
        .when(hl.is_missing(alleles_expr) & ref_is_effect, 2) \
        .when(hl.is_missing(alleles_expr) & ~ref_is_effect, 0) \
        .default(
            hl.or_else(alleles_expr, hl.empty_array(hl.tstr)).filter(
                lambda allele: allele == effect_allele
            ).length()
        )


def _prepare_mt_split(
    vds: hl.vds.VariantDataset,
    weights_table: hl.Table,
    ref_is_effect_allele: bool,
    detailed_timings: bool,
) -> hl.MatrixTable:
    """
    Prepares a MatrixTable for the split-multi PRS calculation path.

    This function takes an interval-filtered VDS, splits multi-allelic sites,
    and joins it with the weights table using a precise `(locus, alleles)`
    key. It handles allele orientation and weight direction based on the
    `ref_is_effect_allele` flag and calculates the bi-allelic dosage.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        The interval-filtered Variant Dataset.
    weights_table : hail.Table
        The chunk of the weights table.
    ref_is_effect_allele : bool
        If True, assumes the effect allele is the reference.
    detailed_timings : bool
        If True, logs the duration of computational steps.

    Returns
    -------
    hail.MatrixTable
        A prepared MatrixTable, filtered and annotated with `weights_info`
        and `dosage` for the specified effect allele.
    """
    with _log_timing(
        "Splitting multi-allelic variants and joining", detailed_timings
    ):
        mt = hl.vds.split_multi(vds).variant_data

        weights_ht_processed = weights_table.annotate(
            alleles=hl.if_else(
                ref_is_effect_allele,
                [weights_table.effect_allele, weights_table.noneffect_allele],
                [weights_table.noneffect_allele, weights_table.effect_allele],
            ),
            weight=hl.if_else(
                ref_is_effect_allele,
                -weights_table.weight,
                weights_table.weight,
            ),
        ).key_by('locus', 'alleles')

        mt = mt.annotate_rows(
            weights_info=weights_ht_processed[mt.row_key]
        )
        mt = mt.filter_rows(hl.is_defined(mt.weights_info))

    with _log_timing(
        "Planning: Calculating per-variant dosage",
        detailed_timings,
    ):
        # After splitting, LGT is converted to GT, so we can
        # directly and safely use the built-in dosage calculator.
        # See the source code for `hl.vds.split_multi` for details.
        mt = mt.annotate_entries(dosage=mt.GT.n_alt_alleles())

        return mt


def _prepare_mt_non_split(
    vds: hl.vds.VariantDataset,
    weights_table: hl.Table,
    strict_allele_match: bool,
    detailed_timings: bool,
) -> hl.MatrixTable:
    """
    Prepares a MatrixTable for the non-split PRS calculation path.

    This function takes an interval-filtered VDS and joins it with the
    weights table using a locus-based key. It optionally performs a strict
    allele match to handle allele orientation and then calculates dosage
    using the custom multi-allelic dosage function.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        The interval-filtered Variant Dataset.
    weights_table : hail.Table
        The chunk of the weights table.
    strict_allele_match : bool
        If True, performs a robust check for allele correspondence.
    detailed_timings : bool
        If True, logs the duration of computational steps.

    Returns
    -------
    hail.MatrixTable
        A prepared MatrixTable, filtered and annotated with `weights_info`
        and `dosage` for the specified effect allele.
    """
    mt = vds.variant_data

    with _log_timing(
        "Planning: Annotating variants with weights", detailed_timings
    ):
        mt = mt.annotate_rows(weights_info=weights_table[mt.locus])
        mt = mt.filter_rows(hl.is_defined(mt.weights_info))

    if strict_allele_match:
        with _log_timing("Performing strict allele match", detailed_timings):
            alt_alleles = hl.set(mt.alleles[1:])
            ref_allele = mt.alleles[0]
            effect = mt.weights_info.effect_allele
            noneffect = mt.weights_info.noneffect_allele

            is_valid_pair = (
                (effect == ref_allele) & alt_alleles.contains(noneffect)
            ) | (
                (noneffect == ref_allele) & alt_alleles.contains(effect)
            )

            mt = mt.filter_rows(is_valid_pair)

    with _log_timing(
        "Planning: Calculating per-variant dosage",
        detailed_timings,
    ):
        mt = mt.annotate_entries(dosage=_calculate_dosage(mt))

    return mt


def _calculate_prs_chunk(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    config: PRSConfig,
) -> hl.Table:
    """
    Calculates a Polygenic Risk Score (PRS) for a single chunk of variants.

    This function serves as the core computational engine. It first filters the
    VDS to a small genomic region based on the input `weights_table` chunk. It
    then dispatches to the appropriate helper function to prepare a MatrixTable
    based on the `split_multi` setting, and finally aggregates the scores.


    Parameters
    ----------
    weights_table : hail.Table
        A pre-filtered chunk of the main weights table, keyed by 'locus'.
    vds : hail.vds.VariantDataset
        The Variant Dataset containing the genetic data.
    config : PRSConfig
        A configuration object containing all settings for the calculation,
        such as `split_multi` and `include_n_matched`.

    Returns
    -------
    hail.Table
        A Hail Table containing the partial PRS results for the chunk, with
        columns for the sample ID, 'prs' score, and optionally 'n_matched'.
    """
    # with _log_timing(
    #     "Planning: Filtering VDS to variants in weights table chunk",
    #     config.detailed_timings,
    # ):
    #     intervals_to_filter = (
    #         weights_table.select(
    #             interval=hl.interval(
    #                 weights_table.locus,
    #                 weights_table.locus,
    #                 includes_end=True,
    #             )
    #         )
    #         .key_by('interval')
    #         .distinct()
    #     )
    #     vds = hl.vds.filter_intervals(vds, intervals_to_filter, keep=True)

    if config.split_multi:
        mt = _prepare_mt_split(
            vds=vds,
            weights_table=weights_table,
            ref_is_effect_allele=config.ref_is_effect_allele,
            detailed_timings=config.detailed_timings,
        )
    else:
        mt = _prepare_mt_non_split(
            vds=vds,
            weights_table=weights_table,
            strict_allele_match=config.strict_allele_match,
            detailed_timings=config.detailed_timings,
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

    prs_table = prs_table.rename({'s': config.sample_id_col})
    # Dropping all global annotations by passing no arguments to select_globals().
    # This is intentional to reduce memory usage.
    return prs_table.select_globals()


def _prepare_weights_for_chunking(
    weights_table: hl.Table,
    weight_col_name: str,
    log_transform_weight: bool,
    chunk_size: typing.Optional[int],
    detailed_timings: bool,
    validate_table: bool = True
) -> tuple[hl.Table, int]:
    """
    Prepares and annotates a weights table for chunked processing.

    This helper function takes the raw weights table, validates it, calculates
    the number of chunks based on the `chunk_size`, and adds a `chunk_id`
    column to each row. This prepares the table for iterative processing in
    the main PRS calculation loop.

    Parameters
    ----------
    weights_table : hail.Table
        The raw input weights table from the user.
    weight_col_name : str
        The name of the column containing effect weights.
    log_transform_weight : bool
        If True, log-transforms the weight column.
    chunk_size : int, optional
        The desired number of variants per chunk. If None, the entire table
        will be processed as a single chunk.
    detailed_timings : bool
        If True, logs the duration of this preparation step.
    validate_table: bool, default True
        If True, the function calls `_validate_and_prepare_weights_table`. If
        False, this validation step is skipped, assuming the table is already
        prepared. This argument is intended for PRS batch calculations, where
        the table validation has already been performed upstream.

    Returns
    -------
    tuple[hail.Table, int]
        A tuple containing:
        - The fully validated and prepared weights table, now annotated with a
          `chunk_id` for each row.
        - An integer representing the total number of chunks.

    Raises
    ------
    ValueError
        If the `weights_table` is empty after the initial validation and
        filtering steps.
    """
    with _log_timing(
        "Preparing and analyzing weights table", detailed_timings
    ):
        if validate_table:
            full_weights_table = _validate_and_prepare_weights_table(
                weights_table, weight_col_name, log_transform_weight
            )
        else:
            full_weights_table = weights_table
        total_variants = full_weights_table.count()
        if total_variants == 0:
            raise ValueError("Weights table is empty after validation.")

        effective_chunk_size = chunk_size or total_variants
        n_chunks = ceil(total_variants / effective_chunk_size)
        logger.info(
            "Total variants: %d, Number of chunks: %d",
            total_variants,
            n_chunks,
        )

        # Don't use chain (hl.Table.add_index().annotate()) as it is not find
        # the idx at the annotation step due to lazy eval.
        full_weights_table = full_weights_table.add_index()
        full_weights_table = full_weights_table.annotate(
            chunk_id=hl.int(
                full_weights_table.idx / effective_chunk_size
            )
        )
        # Note:add_index does not reset the existing key (locus), so don't have
        # to key_by('locus) again

        return full_weights_table, n_chunks


def _process_chunks(
    full_weights_table: hl.Table,
    n_chunks: int,
    vds: hl.vds.VariantDataset,
    config: PRSConfig,
) -> list[pd.DataFrame]:
    """
    Iteratively processes each chunk of the weights table.

    This helper function orchestrates the main processing loop. It iterates
    through the weights table chunk by chunk, calling the core
    `_calculate_prs_chunk` engine for each one, and converts the partial
    results into Pandas DataFrames for final aggregation.

    Parameters
    ----------
    full_weights_table : hail.Table
        The prepared and chunk-annotated weights table.
    n_chunks : int
        The total number of chunks to process.
    vds : hail.vds.VariantDataset
        The Variant Dataset, potentially pre-filtered for samples.
    config : PRSConfig
        The configuration object containing all calculation settings.

    Returns
    -------
    list[pd.DataFrame]
        A list of Pandas DataFrames, where each DataFrame contains the partial
        PRS results for one chunk.
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
            intervals_to_filter = (
                weights_chunk.select(
                    interval=hl.interval(
                        weights_chunk.locus,
                        weights_chunk.locus,
                        includes_end=True,
                    )
                )
                .key_by('interval')
                .distinct()
            )
            vds = hl.vds.filter_intervals(vds, intervals_to_filter, keep=True)

            chunk_prs_table = _calculate_prs_chunk(
                weights_table=weights_chunk,
                vds=vds,
                config=config
            )

            partial_dfs.append(chunk_prs_table.to_pandas())
    return partial_dfs


def _aggregate_and_export(
    partial_dfs: list[pd.DataFrame],
    output_path: str,
    sample_id_col: str,
    detailed_timings: bool,
) -> None:
    """
    Aggregates partial Pandas DataFrame results and exports to a final file.

    This helper function handles the final aggregation and export stage. It
    takes a list of Pandas DataFrames, each containing partial results from a
    single chunk, concatenates them, and then calculates the final total PRS
    for each sample by grouping and summing the results. The final aggregated
    data is then written to a specified cloud storage path.

    Parameters
    ----------
    partial_dfs : list[pd.DataFrame]
        A list of Pandas DataFrames, each containing partial PRS results.
    output_path : str
        The GCS path to write the final tab-separated file.
    sample_id_col : str
        The name of the column containing sample IDs to group by.
    detailed_timings : bool
        If True, logs the duration of the aggregation and export steps.

    Returns
    -------
    None
    """
    if not partial_dfs:
        logger.warning(
            "No PRS results were generated. No output file will be created."
        )
        return

    with _log_timing("Aggregating results with Pandas", detailed_timings):
        combined_df = pd.concat(partial_dfs, ignore_index=True)
        final_df = combined_df.groupby(sample_id_col).sum()

    with _log_timing(
            f"Exporting final result to {output_path}", detailed_timings
    ):
        with hfs.open(output_path, 'w') as f:
            final_df.to_csv(f, sep='\t', index=True, header=True)

def calculate_prs(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    output_path: str,
    config: PRSConfig = PRSConfig(),
) -> typing.Optional[str]:
    """
    Calculates a Polygenic Risk Score (PRS) and exports it to a file.

    This function is the main entry point for the PRS calculation workflow. It
    processes a weights table in chunks, using a filter_intervals approach to
    select variants from the VDS for each chunk. The partial results are then
    converted to Pandas DataFrames and aggregated in memory to produce the
    final score file.

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
        The Variant Dataset containing the genetic data.
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
        The path to the final output file. Returns None if no results are
        generated. The output file is a tab-separated text file with the
        following columns:
        - A sample identifier column (named according to `sample_id_col`).
        - 'prs': The calculated Polygenic Risk Score.
        - 'n_matched' (optional): The number of variants used to calculate
          the score, included if `config.include_n_matched` is True.

    Raises
    ------
    ValueError
        If `output_path` is not a valid GCS path, or if the `weights_table`
        is empty after validation.
    TypeError
        If the `config.samples_to_keep` argument is of an unsupported type.
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
                "Filtering to specified samples", config.detailed_timings
            ):
                samples_ht = _prepare_samples_to_keep(config.samples_to_keep)
                vds = hl.vds.filter_samples(vds, samples_ht)

        full_weights_table, n_chunks = _prepare_weights_for_chunking(
            weights_table=weights_table,
            weight_col_name=config.weight_col_name,
            log_transform_weight=config.log_transform_weight,
            chunk_size=config.chunk_size,
            detailed_timings=config.detailed_timings
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
            sample_id_col=config.sample_id_col,
            detailed_timings=config.detailed_timings
        )

    # Report the total time using the duration captured by the context manager
    logger.info(
        "PRS calculation complete. Total time: %.2f seconds.", timer.duration
    )
    return output_path if partial_dfs else None
