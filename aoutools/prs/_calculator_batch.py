"""PRS batch calculator"""

import typing
import logging
import hail as hl
import hailtop.fs as hfs
import pandas as pd
from aoutools._utils.helpers import SimpleTimer
from ._utils import (
    _log_timing,
)
from ._calculator import (
    _validate_and_prepare_weights_table,
    _prepare_weights_for_chunking,
    _calculate_dosage,
    _prepare_samples_to_keep,
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
        for score_name, table in weights_tables_map.items():
            prepared_table = _validate_and_prepare_weights_table(
                table, config.weight_col_name, config.log_transform_weight
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
                processed_ht = ht.annotate(
                    alleles=hl.if_else(
                        config.ref_is_effect_allele,
                        [ht.effect_allele, ht.noneffect_allele],
                        [ht.noneffect_allele, ht.effect_allele]
                    ),
                    weight=hl.if_else(
                        config.ref_is_effect_allele,
                        -ht.weight,
                        ht.weight
                    )
                ).key_by('locus', 'alleles')
                final_prepared_weights[score_name] = processed_ht
    else:
        final_prepared_weights = prepared_weights

    loci_to_keep = hl.Table.union(*all_loci_tables).key_by('locus').distinct()
    return final_prepared_weights, loci_to_keep


# def _calculate_prs_chunk_batch(
#     vds: hl.vds.VariantDataset,
#     weights_tables_map: dict[str, hl.Table],
#     prepared_weights: dict[str, hl.Table],
#     config: PRSConfig,
# ) -> hl.Table:
#     """
#     Calculates all PRS scores for a single chunk of a VDS.
#     """
#     # Step 1: Set up the MatrixTable and its annotation key based on the path.
#     if config.split_multi:
#         mt = hl.vds.split_multi(vds).variant_data
#         mt_key = mt.row_key
#     else:
#         mt = vds.variant_data
#         mt_key = mt.locus

#     # Step 2: Annotate the rows. This is now common to both paths.
#     annotation_exprs = {
#         f'weights_info_{score_name}': (prepared_weights[score_name][mt_key])
#         for score_name in weights_tables_map
#     }
#     mt = mt.annotate_rows(**annotation_exprs)

#     # Step 3: Build the aggregators for each score.
#     score_aggregators = {}
#     for score_name in weights_tables_map:
#         weights_info = mt[f'weights_info_{score_name}']

#         # The dosage calculation is the main difference between the paths.
#         if config.split_multi:
#             dosage = mt.GT.n_alt_alleles()
#         else:
#             dosage = _calculate_dosage(mt, score_name)

#         partial_score = hl.if_else(
#             hl.is_defined(weights_info),
#             dosage * weights_info.weight,
#             0.0
#         )
#         score_aggregators[score_name] = hl.agg.sum(partial_score)

#     score_aggregators = {}
#     for score_name in weights_tables_map:
#         weights_info = mt[f'weights_info_{score_name}']

#         is_valid_for_score = hl.is_defined(weights_info)

#         if config.split_multi:
#             dosage = mt.GT.n_alt_alleles()
#         else:
#             dosage = _calculate_dosage(mt, score_name)
#             # For the non-split path, conditionally apply the strict allele match
#             if config.strict_allele_match:
#                 alt_alleles = hl.set(mt.alleles[1:])
#                 ref_allele = mt.alleles[0]
#                 effect = weights_info.effect_allele
#                 noneffect = weights_info.noneffect_allele

#                 is_valid_pair = (
#                     (effect == ref_allele) & alt_alleles.contains(noneffect)
#                 ) | (
#                     (noneffect == ref_allele) & alt_alleles.contains(effect)
#                 )

#                 # A score is only valid if the weight is defined AND the
#                 # alleles match
#                 is_valid_for_score = is_valid_for_score & is_valid_pair

#         partial_score = hl.if_else(
#             is_valid_for_score,
#             dosage * weights_info.weight,
#             0.0
#         )
#         score_aggregators[score_name] = hl.agg.sum(partial_score)

#     # Step 4: Run all aggregations in a single pass.
#     return mt.select_cols(**score_aggregators).cols().select_globals()

def _calculate_prs_chunk_batch(
    vds: hl.vds.VariantDataset,
    weights_tables_map: dict[str, hl.Table],
    prepared_weights: dict[str, hl.Table],
    config: PRSConfig,
) -> hl.Table:
    """
    Calculates all PRS scores for a single chunk of a VDS.
    """
    # Step 1: Set up the MatrixTable and its annotation key based on the path.
    if config.split_multi:
        mt = hl.vds.split_multi(vds).variant_data
        mt_key = mt.row_key
    else:
        mt = vds.variant_data
        mt_key = mt.locus

    # Step 2: Annotate the rows. This is now common to both paths.
    annotation_exprs = {
        f'weights_info_{score_name}': prepared_weights[score_name][mt_key]
        for score_name in weights_tables_map
    }
    mt = mt.annotate_rows(**annotation_exprs)

    # Step 3: Build the aggregators for each score.
    score_aggregators = {}
    validity_exprs = {}  # Store validity expressions to reuse for n_matched
    for score_name in weights_tables_map:
        weights_info = mt[f'weights_info_{score_name}']

        is_valid_for_score = hl.is_defined(weights_info)

        if config.split_multi:
            dosage = mt.GT.n_alt_alleles()
        else:
            dosage = _calculate_dosage(mt, score_name)
            if config.strict_allele_match:
                alt_alleles = hl.set(mt.alleles[1:])
                ref_allele = mt.alleles[0]
                effect = weights_info.effect_allele
                noneffect = weights_info.noneffect_allele

                is_valid_pair = (
                    (effect == ref_allele) & alt_alleles.contains(noneffect)
                ) | (
                    (noneffect == ref_allele) & alt_alleles.contains(effect)
                )

                is_valid_for_score = is_valid_for_score & is_valid_pair

        validity_exprs[score_name] = is_valid_for_score

        partial_score = hl.if_else(
            is_valid_for_score,
            dosage * weights_info.weight,
            0.0
        )
        score_aggregators[score_name] = hl.agg.sum(partial_score)

    # Step 4: Run the PRS score aggregations.
    prs_table = mt.select_cols(**score_aggregators).cols().select_globals()

    # Step 5: If requested, calculate n_matched using an efficient single pass.
    if config.include_n_matched:
        with _log_timing(
                "Computing shared variants count", config.detailed_timings
        ):
            # Build a dictionary of count aggregators, one for each score
            n_matched_aggregators = {
                f'n_matched_{score_name}': hl.agg.count_where(
                    validity_exprs[score_name]
                )
                for score_name in weights_tables_map
            }
            # Run all count aggregations in a single pass over the rows
            n_matched_counts = mt.aggregate_rows(
                hl.struct(**n_matched_aggregators)
            )
            # Annotate the results table with the global counts
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

            intervals_to_filter = loci_chunk.select(
                interval=hl.interval(
                    loci_chunk.locus,
                    loci_chunk.locus,
                    includes_end=True
                )
            ).key_by('interval')

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
    sample_id_col: str,
    detailed_timings: bool,
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
            "Aggregating batch results with Pandas", detailed_timings
    ):
        combined_df = pd.concat(partial_dfs, ignore_index=True)
        final_df = combined_df.groupby(sample_id_col).sum()

    with _log_timing(
            f"Exporting final result to {output_path}", detailed_timings
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
            weight_col_name='',
            log_transform_weight=False,
            chunk_size=config.chunk_size,
            detailed_timings=config.detailed_timings,
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
            sample_id_col=config.sample_id_col,
            detailed_timings=config.detailed_timings
        )

    logger.info(
        "Batch PRS calculation complete. Total time: %.2f seconds.",
        timer.duration
    )
    return output_path if partial_dfs else None
