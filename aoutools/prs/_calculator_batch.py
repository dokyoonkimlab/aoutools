"""PRS batch calculator"""

import logging

import hail as hl
import hailtop.fs as hfs
import pandas as pd

from aoutools._utils.helpers import SimpleTimer

from ._calculator_utils import (
    _create_1bp_intervals,
    _entry_contribution,
    _group_weights_by_locus,
    _match_weight_at_locus,
    _orient_weight_and_offset,
    _prepare_samples_to_keep,
    _prepare_weights_for_chunking,
    _split_multi_with_total_dosage,
    _unpersist_quietly,
    _validate_and_prepare_weights_table,
)
from ._config import PRSConfig
from ._utils import _log_timing

logger = logging.getLogger(__name__)


def _prepare_batch_weights_data(
    weights_tables_map: dict[str, hl.Table],
    config: PRSConfig,
) -> tuple[dict, hl.Table]:
    """
    Prepares multiple weights tables for batch PRS calculation.

    This function validates and formats each weights table for scoring. It
    also builds a union of all unique loci across tables, which will later be
    used to filter the Variant Dataset (VDS).

    Parameters
    ----------
    weights_tables_map : dict[str, hl.Table]
        A dictionary mapping score names to Hail tables containing PRS weights.
    config : PRSConfig
        The calculation parameters used to format each weights table -- the
        weight column (``weight_col_name``) and its optional log transform --
        plus ``detailed_timings``.

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
    with _log_timing("Preparing all weights tables", config.detailed_timings):
        for score_name, weights_table in weights_tables_map.items():
            prepared_table = _validate_and_prepare_weights_table(
                weights_table=weights_table, config=config
            )
            prepared_weights[score_name] = prepared_table
            all_loci_tables.append(prepared_table.select())

    if not all_loci_tables:
        return {}, None

    final_prepared_weights = {}
    with _log_timing(
        "Grouping weights tables by locus for the join",
        config.detailed_timings,
    ):
        for score_name, ht in prepared_weights.items():
            final_prepared_weights[score_name] = _group_weights_by_locus(ht)

    loci_to_keep = hl.Table.union(*all_loci_tables).key_by("locus").distinct()
    return final_prepared_weights, loci_to_keep


def _build_row_annotations(
    mt_locus: hl.expr.LocusExpression,
    mt_canonical_alleles: hl.expr.ArrayExpression,
    weights_tables_map: dict[str, hl.Table],
    prepared_weights: dict[str, hl.Table],
) -> dict[str, hl.expr.Expression]:
    """
    Builds a dictionary of row annotations for PRS calculation.

    Each annotation includes:
    - `weights_info_{score}`: A struct containing the weights row matched to
      the current MatrixTable row for the given score.
    - `is_valid_{score}`: A boolean expression indicating whether a valid match
      was found for that score.

    The join is the same shuffle-free one the single-score path uses: match by
    locus (a key-prefix join, `prepared_weights[score][mt_locus]`), then by
    alleles locally against the row's minimal `canonical_alleles`. See
    `_group_weights_by_locus` and `_match_weight_at_locus`.

    Parameters
    ----------
    mt_locus : hl.expr.LocusExpression
        The MatrixTable row's locus (its leading key), joined against each
        locus-grouped weights table.
    mt_canonical_alleles : hl.expr.ArrayExpression
        The row's minimal `[ref, alt]` from `_split_multi_with_total_dosage`,
        matched locally against the weights at that locus.
    weights_tables_map : dict[str, hl.Table]
        A dictionary mapping score names to the original weights tables.
    prepared_weights : dict[str, hl.Table]
        A dictionary of validated, locus-grouped weights tables.

    Returns
    -------
    dict[str, hl.expr.Expression]
        A dictionary mapping annotation names to Hail expressions, to be used as
        row fields in the MatrixTable.
    """
    annotations = {}
    for score_name in weights_tables_map:
        weights_info_expr = _match_weight_at_locus(
            prepared_weights[score_name][mt_locus].variants,
            mt_canonical_alleles,
        )
        annotations[f"weights_info_{score_name}"] = weights_info_expr
        annotations[f"is_valid_{score_name}"] = hl.is_defined(weights_info_expr)
    return annotations


def _build_prs_entry_term(
    mt: hl.MatrixTable,
    score_name: str,
) -> hl.expr.Float64Expression:
    """
    Builds the per-entry aggregation for one score.

    This sums `weight_per_alt_copy * dosage` over the entries of all valid
    variants. A variant is valid if it is matched in the weights table; which
    dosage is counted depends on the row's orientation -- see
    `_entry_contribution` and `_orient_weight_and_offset`.

    This is only the entry term. Homozygous-reference samples have no entry, so
    the row-level `hom_ref_offset` that reaches them is summed separately over
    rows (see `_build_row_offset_expr` and `_calculate_prs_chunk_batch`) and
    added to this term -- which keeps the heavy entry pass to a single scan.

    Unlike the single-score path, batch mode never filters rows -- it masks with
    `hl.if_else(is_valid, ...)` -- so an unmatched weights row must contribute
    nothing here.

    Parameters
    ----------
    mt : hail.MatrixTable
        A MatrixTable containing genotype data and row-level annotations
        produced by `_build_row_annotations`, including `weights_info_{score}`
        and `is_valid_{score}`.
    score_name : str
        A string identifier for the PRS score to compute. Used to look up the
        relevant annotations.

    Returns
    -------
    hl.expr.Float64Expression
        The per-entry aggregation for the given score (without the offset).
    """
    weights_info = mt[f"weights_info_{score_name}"]
    is_valid = mt[f"is_valid_{score_name}"]

    weight_per_alt_copy, hom_ref_offset, ref_is_effect = (
        _orient_weight_and_offset(mt.canonical_alleles[0], weights_info)
    )
    contribution = _entry_contribution(
        mt, weight_per_alt_copy, hom_ref_offset, ref_is_effect
    )

    # Aggregated over entries: hom-ref samples are absent and never visited.
    return hl.agg.sum(hl.if_else(is_valid, contribution, 0.0))


def _build_row_offset_expr(
    rows: hl.Table,
    score_name: str,
) -> hl.expr.Float64Expression:
    """
    The per-row hom-ref offset for one score, masked to matched rows.

    Returns `hom_ref_offset` where the weights row matched, else `0.0`. Summed
    over rows -- a scan that never touches the entry matrix -- this gives the
    constant `2w` that every sample must receive, including the homozygous-
    reference ones that have no entry. See `_orient_weight_and_offset`.

    Takes the **rows table** (`mt.rows()`), not the MatrixTable, so the returned
    expression is sourced from the same table it is aggregated over.

    Batch mode never filters rows, so the offset must be gated on
    `is_valid_{score}`: an unmatched weights row must contribute no offset, or
    every sample is credited for a variant that is not in the callset.
    """
    weights_info = rows[f"weights_info_{score_name}"]
    is_valid = rows[f"is_valid_{score_name}"]
    _, hom_ref_offset, _ = _orient_weight_and_offset(
        rows.canonical_alleles[0], weights_info
    )
    return hl.if_else(is_valid, hom_ref_offset, 0.0)


def _calculate_prs_chunk_batch(
    vds: hl.vds.VariantDataset,
    weights_tables_map: dict[str, hl.Table],
    prepared_weights: dict[str, hl.Table],
    config: PRSConfig,
) -> pd.DataFrame:
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
        The calculation parameters. ``include_n_matched`` controls whether the
        matched-variant counts are also computed here.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with one row per sample and one column per PRS score. If
        `include_n_matched=True`, additional columns for the number of valid
        variants (e.g., `n_matched_score1`) are included. The chunk's cached
        MatrixTable, if it was persisted, is released before returning.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    # Step 1: Get MatrixTable from VDS, splitting multi-allelics
    with _log_timing(
        "Planning: Splitting multi-allelic variants",
        config.detailed_timings,
    ):
        # Also carries each entry's total non-ref count (`n_non_ref`) through
        # the split, which REF-effect rows need; and leaves
        # `filter_changed_loci` at False (raise) on purpose. See
        # `_split_multi_with_total_dosage`.
        mt = _split_multi_with_total_dosage(vds)

    # Step 2: Annotate MatrixTable rows with weights info and validity masks
    with _log_timing(
        "Planning: Calculating and aggregating PRS scores",
        config.detailed_timings,
    ):
        row_annotations = _build_row_annotations(
            mt.locus,
            mt.canonical_alleles,
            weights_tables_map,
            prepared_weights,
        )
        mt = mt.annotate_rows(**row_annotations)

        # Persist when a row reduction (the offsets, or n_matched) runs in
        # addition to the per-sample scores, so split_multi and the join are
        # computed once and the second pass reads the cache. Skipped when
        # effect_allele_is_alt leaves the score as the only pass. See
        # `_calculator.py` for the single-score twin.
        persisted = not config.effect_allele_is_alt or config.include_n_matched
        if persisted:
            mt = mt.persist()

        # Step 3: Reduce the pure *row* quantities -- each score's hom-ref
        # offset, and (if requested) its matched-variant count -- over the rows
        # table, and localize to scalars. Each is a row reduction, a different
        # aggregation scope from the per-sample scores, so Hail computes it in a
        # second pass over the chunk.
        #
        # `effect_allele_is_alt` skips the offset reduction: if every effect
        # allele is the ALT there are no REF-effect variants and every offset is
        # zero. With no `n_matched` the per-sample scores are then the only
        # pass. Orientation stays per row, so a wrong assertion only shifts
        # every score by a constant (rankings preserved), never reordering the
        # cohort. See `_calculator.py` for the single-score twin. The
        # `n_matched` counts are a separate row reduction, so requesting them
        # keeps a second pass (and the persist above) even with this set.
        rows = mt.rows()
        row_agg_exprs = {}
        if not config.effect_allele_is_alt:
            row_agg_exprs.update(
                {
                    f"total_offset_{score_name}": hl.agg.sum(
                        _build_row_offset_expr(rows, score_name)
                    )
                    for score_name in weights_tables_map
                }
            )
        if config.include_n_matched:
            # Batch mode never filters rows (it masks with `is_valid_*`), so the
            # count must be `count_where(is_valid_{score})`, not a plain row
            # count.
            for score_name in weights_tables_map:
                row_agg_exprs[f"n_matched_{score_name}"] = hl.agg.count_where(
                    rows[f"is_valid_{score_name}"]
                )
        row_aggs = (
            rows.aggregate(hl.struct(**row_agg_exprs))
            if row_agg_exprs
            else None
        )

        # Step 4: One entry pass -- each score's entry term, plus its row-level
        # offset scalar unless the offset pass was skipped.
        if config.effect_allele_is_alt:
            score_aggregators = {
                score_name: _build_prs_entry_term(mt, score_name)
                for score_name in weights_tables_map
            }
        else:
            score_aggregators = {
                score_name: _build_prs_entry_term(mt, score_name)
                + row_aggs[f"total_offset_{score_name}"]
                for score_name in weights_tables_map
            }
        prs_table = mt.select_cols(**score_aggregators).cols().select_globals()

        if config.include_n_matched:
            # Per-chunk scalars, identical for every sample; summed across
            # chunks in the pandas aggregation.
            prs_table = prs_table.annotate(
                **{
                    f"n_matched_{score_name}": row_aggs[
                        f"n_matched_{score_name}"
                    ]
                    for score_name in weights_tables_map
                }
            )

        # Rename sample ID column to the user-defined name.
        prs_table = prs_table.rename({"s": config.sample_id_col})

    # Materialize the small per-sample result, then release the chunk's cache.
    # Nothing else unpersists it, so each persisted chunk would stay pinned in
    # executor memory for the whole run; across many chunks that accumulates
    # into memory pressure and spill that erases the persist's benefit.
    result_df = prs_table.to_pandas()
    if persisted:
        _unpersist_quietly(mt)
    return result_df


def _process_chunks_batch(
    # pylint: disable=too-many-arguments
    # pylint: disable=too-many-positional-arguments
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
        The calculation parameters, including ``sample_id_col`` (the sample-ID
        column used in the output).

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
        # Always show chunk processing time to track progress
        with _log_timing(f"Processing chunk {i + 1}/{n_chunks}", True):
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
            # current chunk's loci. Both sides are keyed on locus -- the
            # weights are grouped one row per locus by
            # _group_weights_by_locus -- so a matched locus keeps its whole
            # variants array.
            chunked_prepared_weights = {
                score_name: table.semi_join(loci_chunk)
                for score_name, table in prepared_weights.items()
            }

            # Returns a materialized per-chunk DataFrame (sample ID already
            # renamed), having released its own cached MatrixTable.
            partial_dfs.append(
                _calculate_prs_chunk_batch(
                    vds_chunk,
                    weights_tables_map,
                    chunked_prepared_weights,
                    config,
                )
            )

            # The result is materialized, so the persisted loci chunk can be
            # released too.
            _unpersist_quietly(loci_chunk)
    return partial_dfs


def _aggregate_and_export_batch(
    partial_dfs: list[pd.DataFrame], output_path: str, config: PRSConfig
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
        with hfs.open(output_path, "w") as f:
            final_df.to_csv(f, sep=",", index=True, header=True)


def calculate_prs_batch(
    weights_tables_map: dict[str, hl.Table],
    vds: hl.vds.VariantDataset,
    output_path: str,
    config: PRSConfig | None = None,
) -> str | None:
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
    # PRSConfig is mutable, so a shared default instance would leak state
    # across calls; build a fresh one per call instead.
    config = PRSConfig() if config is None else config

    timer = SimpleTimer()
    with timer:
        if not output_path.startswith("gs://"):
            raise ValueError(
                "The 'output_path' must be a Google Cloud Storage (GCS) "
                "path, starting with 'gs://'."
            )

        logger.info(
            "Starting batch PRS calculation. Final result will be at: %s",
            output_path,
        )

        if config.samples_to_keep is not None:
            with _log_timing(
                "Filtering to specified samples", config.detailed_timings
            ):
                samples_ht = _prepare_samples_to_keep(config.samples_to_keep)
                vds = hl.vds.filter_samples(vds, samples_ht)

        # Step 1: Prepare all weights data and get unique loci
        prepared_weights, loci_to_keep = _prepare_batch_weights_data(
            weights_tables_map, config
        )

        if loci_to_keep is None:
            logger.warning("No variants found in any weights table. Aborting.")
            return None

        # Step 2: Prepare loci for chunked processing
        count = loci_to_keep.count()
        logger.info("Found %d total unique variants across all scores.", count)

        chunked_loci, n_chunks = _prepare_weights_for_chunking(
            weights_table=loci_to_keep, config=config, validate_table=False
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
            partial_dfs=partial_dfs, output_path=output_path, config=config
        )

    logger.info(
        "Batch PRS calculation complete. Total time: %.2f seconds.",
        timer.duration,
    )
    return output_path if partial_dfs else None
