""""PRS batch calculator"""

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
from ._calculator import (
    _validate_and_prepare_weights_table,
    _prepare_mt_split,
    _prepare_mt_non_split,
    _prepare_samples_to_keep,
    _prepare_weights_for_chunking,
    _calculate_dosage,
)
from ._config import PRSConfig

logger = logging.getLogger(__name__)

# NEW HELPER FUNCTION FOR BATCH EXPORT
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

    with _log_timing("Aggregating batch results with Pandas", detailed_timings):
        combined_df = pd.concat(partial_dfs, ignore_index=True)
        final_df = combined_df.groupby(sample_id_col).sum()

    with _log_timing(f"Exporting final result to {output_path}", detailed_timings):
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
    memory-efficient, per-score annotation approach. Supports both
    split-multi and non-split-multi calculation paths.
    """
    timer = SimpleTimer()
    with timer:
        if not output_path.startswith('gs://'):
            raise ValueError(
                "The 'output_path' must be a Google Cloud Storage (GCS) "
                "path, starting with 'gs://'."
            )

        logger.info(
            f"Starting batch PRS calculation (split_multi={config.split_multi}). "
            f"Final result will be at: {output_path}",
        )

        # Prepare each weights table and collect all unique loci
        prepared_weights = {}
        all_loci_tables = []
        with _log_timing("Preparing all weights tables", config.detailed_timings):
            for score_name, table in weights_tables_map.items():
                prepared_table = _validate_and_prepare_weights_table(
                    table, config.weight_col_name, config.log_transform_weight
                )
                prepared_weights[score_name] = prepared_table
                all_loci_tables.append(prepared_table.select())

        # Create joinable tables for the split-multi path if needed
        prepared_weights_split_joinable = {}
        if config.split_multi:
            with _log_timing("Preparing weights tables for split-multi join", config.detailed_timings):
                for score_name, ht in prepared_weights.items():
                    # This logic correctly orients the alleles and weight based on the config flag,
                    # creating a canonical [ref, alt] representation for the join key.
                    processed_ht = ht.annotate(
                        alleles=hl.if_else(
                            config.ref_is_effect_allele,
                            [ht.effect_allele, ht.noneffect_allele],  # [ref, alt]
                            [ht.noneffect_allele, ht.effect_allele]   # [ref, alt]
                        ),
                        # The final_weight now *always* corresponds to the alternate allele.
                        final_weight=hl.if_else(
                            config.ref_is_effect_allele,
                            -ht.weight,  # Flip sign because original weight was for ref
                            ht.weight    # Keep sign because original weight was for alt
                        )
                    ).key_by('locus', 'alleles')
                    prepared_weights_split_joinable[score_name] = processed_ht

        # Get all unique loci and prepare for chunking
        with _log_timing("Finding all unique variants", config.detailed_timings):
            if not all_loci_tables:
                logger.warning("No variants found in any weights table. Aborting.")
                return None

            loci_to_keep = hl.Table.union(*all_loci_tables).key_by('locus').distinct()
            count = loci_to_keep.count()
            logger.info(f"Found {count} total unique variants across all scores.")

            chunked_loci, n_chunks = _prepare_weights_for_chunking(
                weights_table=loci_to_keep,
                weight_col_name='',
                log_transform_weight=False,
                chunk_size=config.chunk_size,
                detailed_timings=config.detailed_timings,
                validate_table=False
            )

        # Process chunks
        partial_dfs = []
        for i in range(n_chunks):
            with _log_timing(f"Processing chunk {i + 1}/{n_chunks}", config.detailed_timings):
                loci_chunk = chunked_loci.filter(chunked_loci.chunk_id == i).persist()
                
                intervals_to_filter = loci_chunk.select(
                    interval=hl.interval(loci_chunk.locus, loci_chunk.locus, includes_end=True)
                ).key_by('interval')
                
                vds_chunk = hl.vds.filter_intervals(vds, intervals_to_filter, keep=True)

                if config.split_multi:
                    # --- SPLIT-MULTI PATH ---
                    split_vds = hl.vds.split_multi(vds_chunk)
                    mt = split_vds.variant_data

                    # Annotate using the pre-keyed (locus, alleles) tables
                    annotation_exprs = {
                        f'weights_info_{score_name}': prepared_weights_split_joinable[score_name][mt.row_key]
                        for score_name in weights_tables_map
                    }
                    mt = mt.annotate_rows(**annotation_exprs)
                    
                    # Build aggregators for each score
                    score_aggregators = {}
                    for score_name in weights_tables_map:
                        weights_info = mt[f'weights_info_{score_name}']
                        
                        # Dosage is simple after splitting
                        dosage = mt.GT.n_alt_alleles()

                        partial_score = hl.if_else(
                            hl.is_defined(weights_info),
                            dosage * weights_info.final_weight, # Use the pre-calculated weight
                            0.0
                        )
                        score_aggregators[score_name] = hl.agg.sum(partial_score)

                else:
                    # --- NON-SPLIT-MULTI PATH ---
                    mt = vds_chunk.variant_data
                    annotation_exprs = {
                        f'weights_info_{score_name}': prepared_weights[score_name][mt.locus]
                        for score_name in weights_tables_map
                    }
                    mt = mt.annotate_rows(**annotation_exprs)

                    score_aggregators = {}
                    for score_name in weights_tables_map:
                        dosage = _calculate_dosage(mt, score_name)
                        partial_score = hl.if_else(
                            hl.is_defined(mt[f'weights_info_{score_name}']),
                            dosage * mt[f'weights_info_{score_name}'].weight,
                            0.0
                        )
                        score_aggregators[score_name] = hl.agg.sum(partial_score)

                # Run all aggregations in a single pass
                chunk_prs_table = mt.select_cols(**score_aggregators).cols()
                chunk_prs_table = chunk_prs_table.rename({'s': config.sample_id_col})
                partial_dfs.append(chunk_prs_table.to_pandas())

        # Aggregate all partial results and write the final file
        _aggregate_and_export_batch(
            partial_dfs=partial_dfs,
            output_path=output_path,
            sample_id_col=config.sample_id_col,
            detailed_timings=config.detailed_timings
        )

    logger.info(
        "Batch PRS calculation complete. Total time: %.2f seconds.", timer.duration
    )
    return output_path if partial_dfs else None
