"""PRS batch calculator"""

from typing import Optional
import logging
import hail as hl
import hailtop.fs as hfs
import pandas as pd
from aoutools._utils.helpers import SimpleTimer
from ._utils import _log_timing
from ._calculator_utils import (
    _prepare_samples_to_keep,
    _validate_and_prepare_weights_table,
    _orient_weights_for_split,
    _check_allele_match,
    _calculate_dosage,
    _prepare_weights_for_chunking,
    _create_1bp_intervals,
)
from ._config import PRSConfig

logger = logging.getLogger(__name__)


def _prepare_batch_weights_data(
    weights_tables_map: dict[str, hl.Table],
    config: PRSConfig,
) -> tuple[dict, hl.Table]:
    """
    Prepares multiple weights tables for batch PRS calculation.

    This function validates and formats each weights table according to the
    selected calculation mode (split or non-split). It also builds a union of
    all unique loci across tables, which will later be used to filter the
    Variant Dataset (VDS).

    Parameters
    ----------
    weights_tables_map : dict[str, hl.Table]
        A dictionary mapping score names to Hail tables containing PRS weights.
    config : PRSConfig
        A configuration object controlling behavior such as whether to split
        multi-allelic variants.

    Returns
    -------

    tuple[dict, hl.Table]
        A tuple containing:
        - A dictionary of prepared weights tables formatted for PRS calculation.
        - A Hail table containing all unique loci to keep for filtering.
          Returns `None` if no tables are provided.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    prepared_weights = {}
    all_loci_tables = []
    with _log_timing(
        "Preparing all weights tables", config.detailed_timings
    ):
        for score_name, weights_table in weights_tables_map.items():
            prepared_table = _validate_and_prepare_weights_table(
                weights_table=weights_table,
                config=config
            )
            prepared_weights[score_name] = prepared_table
            all_loci_tables.append(prepared_table.select())

    if not all_loci_tables:
        return {}, None

    # If split_multi, create a new dictionary with tables formatted for that
    # path. Otherwise, the already prepared tables are used.
    if config.split_multi:
        final_prepared_weights = {}
        with _log_timing(
            "Re-keying weights tables for split-multi join",
            config.detailed_timings
        ):
            for score_name, ht in prepared_weights.items():
                final_prepared_weights[score_name] = \
                    _orient_weights_for_split(ht, config)
    else:
        final_prepared_weights = prepared_weights

    loci_to_keep = hl.Table.union(*all_loci_tables).key_by('locus').distinct()
    return final_prepared_weights, loci_to_keep


def _build_row_annotations(
    mt: hl.MatrixTable,
    mt_key: hl.expr.StructExpression,
    weights_tables_map: dict[str, hl.Table],
    prepared_weights: dict[str, hl.Table],
    config: PRSConfig,
) -> dict[str, hl.expr.Expression]:
    """
    Builds a dictionary of row annotations for PRS calculation.

    Each annotation includes:
    - `weights_info_{score}`: A struct containing the weights row matched to
      the current MatrixTable row for the given score.
    - `is_valid_{score}`: A boolean expression indicating whether a valid match
      was found for that score.

    In non-split mode, if `strict_allele_match` is enabled, the function also
    verifies allele consistency between the MatrixTable and weights table.

    Parameters
    ----------
    mt : hail.MatrixTable
        A MatrixTable containing genotype data.
    mt_key : hl.expr.StructExpression
        A struct expression (e.g., `hl.struct(locus=..., alleles=...)`) used to
        join against each weights table.
    weights_tables_map : dict[str, hl.Table]
        A dictionary mapping score names to the original weights tables.
    prepared_weights : dict[str, hl.Table]
        A dictionary of validated and possibly re-keyed weights tables,
        prepared for lookup during annotation.
    config : PRSConfig
        A configuration object that controls behavior such as whether to split
        multi-allelic variants and enforce strict allele matching.

    Returns
    -------
    dict[str, hl.expr.Expression]
        A dictionary mapping annotation names to Hail expressions, to be used as
        row fields in the MatrixTable.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    annotations = {}
    for score_name in weights_tables_map:
        weights_info_expr = prepared_weights[score_name][mt_key]
        is_valid_expr = hl.is_defined(weights_info_expr)

        if not config.split_multi and config.strict_allele_match:
            is_valid_expr &= _check_allele_match(mt, weights_info_expr)

        annotations[f'weights_info_{score_name}'] = weights_info_expr
        annotations[f'is_valid_{score_name}'] = is_valid_expr
    return annotations


def _build_prs_agg_expr(
    mt: hl.MatrixTable,
    score_name: str,
    config: PRSConfig,
) -> hl.expr.Aggregation:
    """
    Builds an aggregation expression to compute a Polygenic Risk Score (PRS)
    for a given score.

    This expression calculates the sum of dosage multiplied by weight across
    all valid variants. A variant is considered valid if it is matched in the
    weights table and, if applicable, passes the allele-matching check.

    The dosage computation differs based on whether multi-allelic variants have
    been split. In split mode, `n_alt_alleles()` is used directly. In non-split
    mode, dosage is computed using the custom logic defined in
    `_calculate_dosage`.

    Parameters
    ----------
    mt : hail.MatrixTable
        A MatrixTable containing genotype data and row-level annotations
        produced by `_build_row_annotations`, including `weights_info_{score}`
        and `is_valid_{score}`.
    score_name : str
        A string identifier for the PRS score to compute. Used to look up the
        relevant annotations.
    config : PRSConfig
        A configuration object that determines whether the input data is split
        and how dosage is calculated.

    Returns
    -------
    hl.expr.Aggregation
        An aggregation expression representing the sum of `dosage * weight`
        across all valid variants for the given score.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    weights_info = mt[f'weights_info_{score_name}']
    is_valid = mt[f'is_valid_{score_name}']

    # Use appropriate dosage calculation based on whether VDS is split
    if config.split_multi:
        dosage = mt.GT.n_alt_alleles()
    else:
        dosage = _calculate_dosage(mt, score_name)

    # Final score = sum of (dosage * weight) for valid variants
    return hl.agg.sum(hl.if_else(is_valid, dosage * weights_info.weight, 0.0))


def _calculate_prs_chunk_batch(
    vds: hl.vds.VariantDataset,
    weights_tables_map: dict[str, hl.Table],
    prepared_weights: dict[str, hl.Table],
    config: PRSConfig,
) -> hl.Table:
    """
    Calculates all Polygenic Risk Scores (PRS) for a single chunk of a
    VariantDataset (VDS).

    This function processes a subset of variants from a VDS and computes PRS
    values for all configured scores. It handles variant splitting, annotation
    with weight information, and dosage-based aggregation per sample.

    If configured, it also computes the number of valid variants (i.e., matched
    between the VDS and each weights table) used in score calculation.

    Parameters
    ----------
    vds : hl.vds.VariantDataset
        A VariantDataset chunk containing genotype and variant information.
    weights_tables_map : dict[str, hl.Table]
        A dictionary mapping score names to their original weights tables.
    prepared_weights : dict[str, hl.Table]
        A dictionary of weights tables that have been validated and formatted
        for PRS computation.
    config : PRSConfig
        A configuration object that controls behavior such as whether to split
        multi-allelic variants, whether to perform strict allele matching, and
        whether to include matched variant counts.

    Returns
    -------
    hl.Table
        A Hail Table with one row per sample and one column per PRS score. If
        `include_n_matched=True`, additional columns for the number of valid
        variants (e.g., `n_matched_score1`) are included.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    # Step 1: Get MatrixTable from VDS, optionally splitting multi-allelics
    if config.split_multi:
        with _log_timing(
                "Planning: Splitting multi-allelic variants",
                config.detailed_timings
        ):
            mt = hl.vds.split_multi(vds).variant_data
            mt_key = mt.row_key
    else:
        mt = vds.variant_data
        mt_key = mt.locus

    # Step 2: Annotate MatrixTable rows with weights info and validity masks
    with _log_timing(
        "Planning: Calculating and aggregating PRS scores",
        config.detailed_timings
    ):
        row_annotations = _build_row_annotations(
            mt, mt_key, weights_tables_map, prepared_weights, config
        )
        mt = mt.annotate_rows(**row_annotations)

        # Step 3: Build score aggregators across columns (i.e., samples)
        score_aggregators = {
            score_name: _build_prs_agg_expr(mt, score_name, config)
            for score_name in weights_tables_map
        }

        # Compute and return the per-sample PRS results
        prs_table = mt.select_cols(**score_aggregators).cols().select_globals()

    # Step 4 (Optional): Compute number of matched variants (n_matched_*) if
    # requested
    if config.include_n_matched:
        with _log_timing(
                "Computing shared variants count", config.detailed_timings
        ):
            n_matched_aggs = {
                f'n_matched_{score_name}': hl.agg.count_where(
                    mt[f'is_valid_{score_name}']
                )
                for score_name in weights_tables_map
            }
            n_matched_counts = mt.aggregate_rows(hl.struct(**n_matched_aggs))
            prs_table = prs_table.annotate(**n_matched_counts)

    return prs_table


def _process_chunks_batch(
    #pylint: disable=too-many-arguments
    #pylint: disable=too-many-positional-arguments
    n_chunks: int,
    chunked_loci: hl.Table,
    vds: hl.vds.VariantDataset,
    weights_tables_map: dict[str, hl.Table],
    prepared_weights: dict[str, hl.Table],
    config: PRSConfig,
) -> list[pd.DataFrame]:
    """
    Processes each genomic chunk to compute PRS values in batch mode.

    This function iterates over genomic chunks defined in `chunked_loci`,
    filters the VariantDataset (VDS) to each chunk, and calculates Polygenic
    Risk Scores (PRS) for all configured scores. The results from each chunk
    are returned as a list of pandas DataFrames, one per chunk.

    Parameters
    ----------
    n_chunks : int
        The total number of genomic chunks to process.
    chunked_loci : hl.Table
        A Hail Table containing loci grouped by chunk ID.
    vds : hl.vds.VariantDataset
        A VariantDataset containing the full set of genotypes.
    weights_tables_map : dict[str, hl.Table]
        A dictionary mapping score names to their original weights tables.
    prepared_weights : dict[str, hl.Table]
        A dictionary of weights tables that have been validated and formatted
        for PRS computation.
    config : PRSConfig
        A configuration object that controls behavior such as whether to split
        multi-allelic variants and which sample ID column to use.

    Returns
    -------
    list[pd.DataFrame]
        A list of pandas DataFrames, one per chunk, each containing per-sample
        PRS results and optionally matched variant counts.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    partial_dfs = []
    for i in range(n_chunks):
        with _log_timing(
            f"Processing chunk {i + 1}/{n_chunks}",
            config.detailed_timings
        ):
            loci_chunk = chunked_loci.filter(
                chunked_loci.chunk_id == i
            ).persist()

            intervals_to_filter = _create_1bp_intervals(loci_chunk)
            # If filter_intervals filters the main vds and reassigns to vds
            # again, subsequent operation will try to filter empty variable.
            vds_chunk = hl.vds.filter_intervals(
                vds, intervals_to_filter, keep=True
            )

            # Filter the full weights table via semi-join to retain only the
            # current chunk's weights. semi_join preserves split variants at
            # the same locus when split_multi=True.
            chunked_prepared_weights = {
                score_name: table.semi_join(loci_chunk)
                for score_name, table in prepared_weights.items()
            }

            chunk_prs_table = _calculate_prs_chunk_batch(
                vds_chunk,
                weights_tables_map,
                chunked_prepared_weights,
                config
            )

            chunk_prs_table = chunk_prs_table.rename(
                {'s': config.sample_id_col}
            )
            partial_dfs.append(chunk_prs_table.to_pandas())
    return partial_dfs


def _aggregate_and_export_batch(
    partial_dfs: list[pd.DataFrame],
    output_path: str,
    config: PRSConfig
) -> None:
    """
    Aggregates partial PRS results from all chunks and exports to disk.

    This function combines the list of pandas DataFrames produced by
    `_process_chunks_batch`, sums PRS scores across all chunks for each sample,
    and writes the final per-sample results to a comma-delimited file.

    Parameters
    ----------
    partial_dfs : list[pd.DataFrame]
        A list of pandas DataFrames containing chunk-wise PRS results.
    output_path : str
        A destination path on GCS to write the final comma-separated file.
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
            "Aggregating batch results with Pandas", config.detailed_timings
    ):
        combined_df = pd.concat(partial_dfs, ignore_index=True)
        final_df = combined_df.groupby(config.sample_id_col).sum()

    with _log_timing(
            f"Exporting final result to {output_path}", config.detailed_timings
    ):
        with hfs.open(output_path, 'w') as f:
            final_df.to_csv(f, sep=',', index=True, header=True)


def calculate_prs_batch(
    weights_tables_map: dict[str, hl.Table],
    vds: hl.vds.VariantDataset,
    output_path: str,
    config: PRSConfig = PRSConfig(),
) -> Optional[str]:
    """
    Calculates multiple Polygenic Risk Scores (PRS) concurrently using a
    memory-efficient, per-score annotation approach.

    This function performs a batch PRS calculation on a Hail VariantDataset,
    using chunked aggregation and optional sample filtering.

    Parameters
    ----------
    weights_tables_map : dict[str, hl.Table]
        A dictionary mapping score names to their corresponding PRS weights
        tables.
    vds : hl.vds.VariantDataset
        A Hail VariantDataset containing both variant and sample data.
    output_path : str
        A GCS path (starting with 'gs://') to write the final comma-separated
        output file.
    config : PRSConfig, optional
        A configuration object for all optional parameters. If not provided,
        default settings will be used. See the `PRSConfig` class for details
        on all available settings.

    Returns
    -------
    Optional[str]
        The path to the final PRS result file if successful; otherwise, `None`
        if no valid variants were found.

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
                "Filtering to specified samples", config.detailed_timings
            ):
                samples_ht = _prepare_samples_to_keep(config.samples_to_keep)
                vds = hl.vds.filter_samples(vds, samples_ht)

        # Step 1: Prepare all weights data and get unique loci
        prepared_weights, loci_to_keep = \
            _prepare_batch_weights_data(weights_tables_map, config)

        if loci_to_keep is None:
            logger.warning(
                "No variants found in any weights table. Aborting."
            )
            return None

        # Step 2: Prepare loci for chunked processing
        count = loci_to_keep.count()
        logger.info(
            "Found %d total unique variants across all scores.", count
        )

        chunked_loci, n_chunks = _prepare_weights_for_chunking(
            weights_table=loci_to_keep,
            config=config,
            validate_table=False
        )

        # Step 3: Process all chunks using the new helper function
        partial_dfs = _process_chunks_batch(
            n_chunks=n_chunks,
            chunked_loci=chunked_loci,
            vds=vds,
            weights_tables_map=weights_tables_map,
            prepared_weights=prepared_weights,
            config=config,
        )

        # Step 4: Aggregate and export final results
        _aggregate_and_export_batch(
            partial_dfs=partial_dfs,
            output_path=output_path,
            config=config
        )

    logger.info(
        "Batch PRS calculation complete. Total time: %.2f seconds.",
        timer.duration
    )
    return output_path if partial_dfs else None
