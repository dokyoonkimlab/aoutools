"""PRS calculator"""

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
)
from ._config import PRSConfig
from ._utils import _log_timing

logger = logging.getLogger(__name__)


def _prepare_mt_split(
    vds: hl.vds.VariantDataset, weights_table: hl.Table, config: PRSConfig
) -> hl.MatrixTable:
    """
    Prepares a MatrixTable for PRS calculation.

    Splits multi-allelic sites in the VDS, matches the weights by locus and
    then by minimal alleles (a shuffle-free join; see `_group_weights_by_locus`
    and `_match_weight_at_locus`), and resolves each matched row against the
    VDS's own REF/ALT orientation. Dosage is `GT.n_alt_alleles()` after
    splitting.

    Rows are annotated with `weight_per_alt_copy` and `hom_ref_offset`. The
    caller must add the offset -- a row-level constant -- to every sample; see
    `_orient_weight_and_offset` for why it cannot be an entry-level term.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        An interval-filtered VariantDataset.
    weights_table : hail.Table
        A chunk of the PRS weights table.
    config : PRSConfig
        A configuration object controlling PRS behavior, including
        `detailed_timings`.

    Returns
    -------
    hail.MatrixTable
        A MatrixTable annotated with `weights_info`, `weight_per_alt_copy`,
        `hom_ref_offset`, and per-variant `dosage`.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS
    calculation.
    """
    with _log_timing(
        "Planning: Splitting multi-allelic variants and joining",
        config.detailed_timings,
    ):
        # Splits multi-allelic sites, and carries each entry's *total* non-ref
        # allele count (`n_non_ref`) through the split. Splitting downcodes
        # other ALTs to the reference, which would otherwise make a sample
        # carrying a different ALT indistinguishable from a hom-ref sample --
        # wrong for any REF-effect variant. See `_split_multi_with_total_dosage`
        # for the mechanism and the deliberate `filter_changed_loci` tripwire.
        mt = _split_multi_with_total_dosage(vds)

        # Match weights by locus (a shuffle-free key-prefix join) and then by
        # alleles locally against the row's minimal `canonical_alleles`. Keying
        # the weights on the full (locus, alleles) instead would force the MT to
        # be re-keyed to minimal alleles -- a shuffle of every entry. See
        # `_group_weights_by_locus` and `_match_weight_at_locus`.
        weights_by_locus = _group_weights_by_locus(weights_table)
        mt = mt.annotate_rows(
            weights_info=_match_weight_at_locus(
                weights_by_locus[mt.locus].variants, mt.canonical_alleles
            )
        )

        # Only rows that matched a variant in the VDS survive. This filter is
        # what makes it safe to add `hom_ref_offset` to every sample below: an
        # unmatched weights row contributes no offset, because its row is gone.
        mt = mt.filter_rows(hl.is_defined(mt.weights_info))

    with _log_timing(
        "Planning: Calculating per-variant dosage",
        config.detailed_timings,
    ):
        weight_per_alt_copy, hom_ref_offset, ref_is_effect = (
            _orient_weight_and_offset(mt.canonical_alleles[0], mt.weights_info)
        )
        mt = mt.annotate_rows(
            weight_per_alt_copy=weight_per_alt_copy,
            hom_ref_offset=hom_ref_offset,
            ref_is_effect=ref_is_effect,
        )

        mt = mt.annotate_entries(
            contribution=_entry_contribution(
                mt,
                mt.weight_per_alt_copy,
                mt.hom_ref_offset,
                mt.ref_is_effect,
            )
        )

        return mt


def _calculate_prs_chunk(
    weights_table: hl.Table, vds: hl.vds.VariantDataset, config: PRSConfig
) -> pd.DataFrame:
    """
    Calculates a Polygenic Risk Score (PRS) for a single chunk of variants.

    This function serves as the core computation step. It splits multi-allelic
    sites, joins the variant data to the weights, and computes the PRS using
    dosage-weight aggregation.

    Parameters
    ----------
    weights_table : hail.Table
        A pre-filtered chunk of the full weights table, keyed by 'locus'.
    vds : hail.vds.VariantDataset
        The Variant Dataset containing genotypes to score.
    config : PRSConfig
        A configuration object specifying settings such as `include_n_matched`
        and `sample_id_col`.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with one row per sample and a PRS column. If requested,
        also includes the number of matched variants ('n_matched'). The
        chunk's cached MatrixTable, if it was persisted, is released before
        returning.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    mt = _prepare_mt_split(
        vds=vds,
        weights_table=weights_table,
        config=config,
    )

    # When a row reduction (the offset, or n_matched) runs in addition to the
    # per-sample score, Hail reads the chunk twice. Persist so split_multi and
    # the join are computed once and the second pass reads the cache -- ~10-20%
    # on a many-chunk run such as a 1M-variant score. With effect_allele_is_alt
    # and no n_matched the score is the only pass, so persisting would just add
    # a wasted materialization.
    persisted = not config.effect_allele_is_alt or config.include_n_matched
    if persisted:
        mt = mt.persist()

    # The hom-ref offset -- `2w` summed over the matched rows whose effect
    # allele is the REF base -- reaches homozygous-reference samples, who have
    # no entry in `variant_data` and so cannot be reached by the per-sample
    # entry aggregation. It is a *row* reduction, a different aggregation scope
    # from the per-sample score, so Hail computes it in a second pass.
    #
    # `effect_allele_is_alt` lets the caller skip that pass. If every variant's
    # effect allele is the ALT (non-reference) base -- the convention of most
    # harmonized GWAS summaries -- then no variant is REF-effect, the offset is
    # identically zero, and its reduction is skipped. With no `n_matched` the
    # per-sample score is then the only pass. Orientation is still resolved per
    # row, so if the assertion is wrong the only effect is that every score is
    # short by the same constant (`2w` per REF-effect variant): absolute values
    # shift, but rankings, percentiles, and z-scores are exact, and the cohort
    # is never reordered. `n_matched` is a separate row reduction, so requesting
    # it keeps a second pass (and the persist above) even with this set.
    total_offset = 0.0
    n_matched = None
    if not config.effect_allele_is_alt:
        rows = mt.rows()
        row_agg_exprs = {"total_offset": hl.agg.sum(rows.hom_ref_offset)}
        if config.include_n_matched:
            # `mt` is filtered to matched rows (`filter_rows` in
            # `_prepare_mt_split`), so a plain row count *is* the number of
            # variants in common. A genotype-aware `hl.agg.count()` over entries
            # would undercount, since hom-ref samples have no entry.
            row_agg_exprs["n_matched"] = hl.agg.count()
        row_aggs = rows.aggregate(hl.struct(**row_agg_exprs))
        total_offset = row_aggs.total_offset
        if config.include_n_matched:
            n_matched = row_aggs.n_matched
    elif config.include_n_matched:
        n_matched = mt.rows().aggregate(hl.agg.count())

    # One pass over entries: the per-sample entry sum plus the row-level offset
    # constant (a scalar, identical for every sample; zero when the offset pass
    # was skipped).
    prs_table = mt.select_cols(
        prs=hl.agg.sum(mt.contribution) + total_offset
    ).cols()

    if config.include_n_matched:
        # A per-chunk scalar; summed across chunks in the pandas aggregation.
        prs_table = prs_table.annotate(n_matched=n_matched)

    # Rename sample ID column to user-defined name
    prs_table = prs_table.rename({"s": config.sample_id_col})

    # Drop all global annotations to minimize memory footprint
    prs_table = prs_table.select_globals()

    # Materialize the small per-sample result, then release the chunk's cache.
    # Nothing else unpersists it, so each persisted chunk would stay pinned in
    # executor memory for the whole run; across many chunks that accumulates
    # into memory pressure and spill that erases the persist's benefit.
    result_df = prs_table.to_pandas()
    if persisted:
        _unpersist_quietly(mt)
    return result_df


def _process_chunks(
    full_weights_table: hl.Table,
    n_chunks: int,
    vds: hl.vds.VariantDataset,
    config: PRSConfig,
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
        # Always show chunk processing time to track progress
        with _log_timing(f"Processing chunk {i + 1}/{n_chunks}", True):
            # Use .persist() to avoid recomputation of the same chunk in
            # _calculate_prs_chunk, specifically during:
            # 1. Creation of interval_ht
            # 2. Annotating rows with PRS weight information
            weights_chunk = full_weights_table.filter(
                full_weights_table.chunk_id == i
            ).persist()

            intervals_to_filter = _create_1bp_intervals(weights_chunk)
            # If filter_intervals filters the main vds and reassigns to vds
            # again, subsequent operation will try to filter empty variable.
            vds_chunk = hl.vds.filter_intervals(
                vds, intervals_to_filter, keep=True
            )

            # Returns a materialized per-chunk DataFrame, having released its
            # own cached MatrixTable.
            partial_dfs.append(
                _calculate_prs_chunk(
                    weights_table=weights_chunk, vds=vds_chunk, config=config
                )
            )

            # The result is materialized, so the persisted weights chunk can
            # be released too.
            _unpersist_quietly(weights_chunk)

    return partial_dfs


def _aggregate_and_export(
    partial_dfs: list[pd.DataFrame], output_path: str, config: PRSConfig
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
        "Aggregating results with Pandas", config.detailed_timings
    ):
        combined_df = pd.concat(partial_dfs, ignore_index=True)
        final_df = combined_df.groupby(config.sample_id_col).sum()

    with _log_timing(
        f"Exporting final result to {output_path}", config.detailed_timings
    ):
        with hfs.open(output_path, "w") as f:
            final_df.to_csv(f, sep=",", index=True, header=True)


def calculate_prs(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    output_path: str,
    config: PRSConfig | None = None,
) -> str | None:
    """
    Calculates a Polygenic Risk Score (PRS) and exports the result to a file.

    This function is the main entry point for the PRS calculation workflow. It
    processes a weights table in chunks, using a filter_intervals approach to
    select variants from the VDS for each chunk. Partial results are then
    converted to Pandas DataFrames and aggregated to produce the final score
    file.

    Notes
    -----
    Multi-allelic variants in the VDS are always split before scoring. Splitting
    puts each allele in its minimal representation, which lets a simple variant
    in the weights table match a complex one in the VDS: a weights row for
    chr1:10075251 A/G matches VDS alleles ['AGGGC', 'A', 'GGGGC'], because
    'AGGGC' -> 'GGGGC' minimizes to ['A', 'G'] at that same locus.

    Minimal representation can also *move* a variant's locus, when the two
    alleles share a leading base: ['GG', 'G', 'GT'] at chr1:1001 minimizes its
    SNP allele to ['G', 'T'] at chr1:**1002**. Hail will not relocate a row, so
    such a variant **raises** rather than being scored. That is deliberate. No
    variant of this shape exists in All of Us -- 0 of 6,001,424 ALT alleles
    measured -- so the exception is a tripwire against a future change in
    variant representation, not a crash risk. See `TODO.md`, Finding 5.

    Effect alleles may be given in either orientation, and are resolved per row
    against the VDS's own REF/ALT. The PGS Catalog does not harmonize the effect
    allele onto the ALT, so a file mixing both is the normal case.

    Parameters
    ----------
    weights_table : hail.Table
        A Hail table containing variant weights. Must contain the following
        columns:

        - `chr`: str
        - `pos`: int32
        - `effect_allele`: str
        - `noneffect_allele`: str
        - A column for the effect weight (float64), specified by
          `weight_col_name`.
    vds : hail.vds.VariantDataset
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
    str or None
        The output path if results are successfully written; otherwise, None.
        The output file is a comma-separated text file with:

        - A sample ID column (as configured in `config.sample_id_col`)
        - `prs`: The calculated PRS value
        - `n_matched` (optional): The number of variants used to calculate
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
            "Starting PRS calculation. Final result will be at: %s",
            output_path,
        )

        if config.samples_to_keep is not None:
            with _log_timing(
                "Planning: Filtering to specified samples",
                config.detailed_timings,
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
