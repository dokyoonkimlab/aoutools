"""PRS batch calculator"""

import typing
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
    Prepares all weights tables for batch processing.

    This function validates each weights table, prepares a version formatted
    for the specified calculation path (split or non-split), and returns the
    union of all unique loci found across all tables.

    Returns
    -------
    tuple[dict, hl.Table]
        A tuple containing:
        - A dictionary of prepared weights, formatted for the chosen path.
        - A Hail Table with all unique loci to keep.
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
                # processed_ht = ht.annotate(
                #     alleles=hl.if_else(
                #         config.ref_is_effect_allele,
                #         [ht.effect_allele, ht.noneffect_allele],
                #         [ht.noneffect_allele, ht.effect_allele]
                #     ),
                #     weight=hl.if_else(
                #         config.ref_is_effect_allele,
                #         -ht.weight,
                #         ht.weight
                #     )
                # ).key_by('locus', 'alleles')
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
    Returns a dictionary of row annotations including:
    - weights_info_{score}
    - is_valid_{score}
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
    Returns an aggregation expression for a given score name.
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
    Calculates all PRS scores for a single chunk of a VDS.
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
    n_chunks: int,
    chunked_loci: hl.Table,
    vds: hl.vds.VariantDataset,
    weights_tables_map: dict[str, hl.Table],
    prepared_weights: dict[str, hl.Table],
    config: PRSConfig,
) -> list[pd.DataFrame]:
    """
    Iteratively processes each chunk of loci for batch PRS calculation.
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
            vds = hl.vds.filter_intervals(vds, intervals_to_filter, keep=True)

            chunk_prs_table = _calculate_prs_chunk_batch(
                vds,
                weights_tables_map,
                prepared_weights,
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
    # sample_id_col: str,
    # detailed_timings: bool,
    config: PRSConfig
) -> None:
    """
    Aggregates partial results from per-score aggregation and exports.
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
            final_df.to_csv(f, sep='\t', index=True, header=True)


def calculate_prs_batch(
    weights_tables_map: dict[str, hl.Table],
    vds: hl.vds.VariantDataset,
    output_path: str,
    config: PRSConfig = PRSConfig(),
) -> typing.Optional[str]:
    """
    Calculates multiple Polygenic Risk Scores (PRS) concurrently using a
    memory-efficient, per-score annotation approach.
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
