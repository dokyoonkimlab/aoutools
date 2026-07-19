"""Internal utilities for the PRS calculator"""

# See the note in _config.py: mocked hail (docs builds) cannot evaluate the
# PEP 604 `|` in these annotations, so defer them.
from __future__ import annotations

import logging
from math import ceil

import hail as hl
from hail.utils.java import FatalError

from ._config import PRSConfig
from ._utils import (
    _log_timing,
    _standardize_chromosome_column,
)

logger = logging.getLogger(__name__)


def _unpersist_quietly(dataset: hl.MatrixTable | hl.Table) -> None:
    """
    Release a persisted Hail dataset's cache, best-effort.

    A persisted chunk that is never released stays pinned in executor memory
    for the whole run, so freeing it once its result is materialized keeps
    memory from growing across chunks. Hail's local Spark backend, however,
    removes the persist directory non-recursively and raises
    ``IOException: Directory ... is not empty`` -- a backend quirk, not a data
    problem. The result is already computed by the time this is called, and any
    cache the release leaves behind is reclaimed when the Hail session ends, so
    a failure here is safe to ignore.
    """
    try:
        dataset.unpersist()
    except FatalError as exc:
        logger.debug("Could not release a persisted cache: %s", exc)


def _prepare_samples_to_keep(
    samples: hl.Table | list | set | tuple | int | float | str,
) -> hl.Table:
    """
    Converts a flexible list of samples into a keyed Hail Table.

    Accepts various common Python collection types (list, set, tuple) or single
    values (int, float, str) and converts them into a standardized Hail Table
    keyed by the string column 's'. This table can then be used to filter
    Hail datasets by sample ID.

    Parameters
    ----------
    samples : hail.Table, list, set, tuple, int, float, or str
        A collection or single value of sample IDs to prepare. Numeric types
        will be converted to strings for consistent keying.

    Returns
    -------
    hail.Table
        A Hail Table keyed by 's' containing sample IDs as strings.

    Raises
    ------
    TypeError
        If the input `samples` object is not one of the supported types.
    """
    if isinstance(samples, hl.Table):
        return samples

    if isinstance(samples, (int, float, str)):
        sample_list = [str(samples)]
    elif isinstance(samples, (list, set, tuple)):
        sample_list = [str(s) for s in samples]
    else:
        raise TypeError(
            f"Unsupported type for samples_to_keep: {type(samples)}."
        )

    samples_ht = hl.Table.parallelize(
        [{"s": s} for s in sample_list], hl.tstruct(s=hl.tstr)
    )
    return samples_ht.key_by("s")


def _validate_and_prepare_weights_table(
    weights_table: hl.Table, config: PRSConfig
) -> hl.Table:
    """
    Validates and prepares a single weights table for PRS calculation.

    Ensures the table has the required columns with correct types, standardizes
    the chromosome format, handles the weight column (renaming and optional
    log transformation), and keys the table by locus for joining with the VDS.

    Parameters
    ----------
    weights_table : hail.Table
        An input weights table. Must contain columns 'chr', 'pos',
        'effect_allele', 'noneffect_allele', and a weight column.
    config : PRSConfig
        A configuration object specifying parameters such as `weight_col_name`
        and `log_transform_weight`.

    Returns
    -------
    hail.Table
        A validated, standardized, and keyed Hail Table.

    Raises
    ------
    TypeError
        If `weight_col_name` is missing, or any required column is missing, or
        has incorrect data type.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    if config.weight_col_name not in weights_table.row:
        raise TypeError(
            f"Specified weight column '{config.weight_col_name}' not found "
            f"in table."
        )
    weights_table = weights_table.rename({config.weight_col_name: "weight"})

    required_cols = {
        "chr": hl.tstr,
        "pos": hl.tint32,
        "effect_allele": hl.tstr,
        "noneffect_allele": hl.tstr,
        "weight": hl.tfloat64,
    }
    for col, expected_type in required_cols.items():
        if col not in weights_table.row:
            raise TypeError(
                f"Weights table is missing required column: '{col}'."
            )
        if weights_table[col].dtype != expected_type:
            raise TypeError(f"Column '{col}' has incorrect type.")

    weights_table = _standardize_chromosome_column(weights_table)

    if config.log_transform_weight:
        weights_table = weights_table.annotate(
            weight=hl.log(weights_table.weight)
        )

    weights_table = weights_table.annotate(
        locus=hl.locus(
            weights_table.chr, weights_table.pos, reference_genome="GRCh38"
        )
    )
    weights_table = weights_table.key_by("locus")
    return weights_table.select("effect_allele", "noneffect_allele", "weight")


def _group_weights_by_locus(ht: hl.Table) -> hl.Table:
    """
    Groups a weights table by locus into a per-locus array of variants.

    Returns a table keyed by `locus` alone, with a `variants` field holding an
    array of every weights row at that locus (`effect_allele`,
    `noneffect_allele`, `weight`). `_match_weight_at_locus` picks the right one.

    Why key by locus and not by `(locus, alleles)`. The split `variant_data`
    is keyed by `(locus, alleles)`, partitioned by locus, but its alleles are
    the VDS's own -- possibly non-minimal -- representation, while a GWAS names
    variants minimally. Matching on the row key therefore requires *rewriting*
    the VDS alleles to minimal form, which means re-keying the MatrixTable, and
    a re-key shuffles the whole chunk (entries and all).

    Keying the weights by locus alone avoids that. `weights_by_locus[mt.locus]`
    is a key-*prefix* join against the MT's leading key, which hail does with no
    shuffle -- only the small weights table is grouped. Alleles are then matched
    locally against `mt.canonical_alleles` (the minimal form annotated by
    `_split_multi_with_total_dosage`), by `_match_weight_at_locus`. So the big
    entry data never moves. This is the whole reason the split MT is left on its
    own key rather than normalized in place.

    The locus alone is not a *variant* identity -- a file may name more than one
    variant at a position -- which is exactly why `variants` is an array and the
    allele match is kept: a weights row for an A/G SNP must not match an A/T
    variant at the same locus. A degenerate row naming the same allele twice
    (`A`/`A`) can match no VDS row (REF != ALT) and is dropped here.

    Orientation is read per row off the VDS's REF (`_orient_weight_and_offset`),
    not declared file-wide, because the PGS Catalog does not harmonize the
    effect allele onto the ALT: its spec says the effect allele "does not
    necessarily need to correspond to the minor allele/alternative allele", and
    harmonized (HmPOS) files remap coordinates while preserving the original
    allele columns. A single file routinely carries both orientations, which the
    unordered allele-set match in `_match_weight_at_locus` handles directly.

    Parameters
    ----------
    ht : hail.Table
        A Hail table containing 'effect_allele', 'noneffect_allele', 'weight',
        and 'locus' fields.

    Returns
    -------
    hail.Table
        A table keyed by 'locus', with a 'variants' array of
        `struct(effect_allele, noneffect_allele, weight)`.
    """
    ht = ht.filter(ht.effect_allele != ht.noneffect_allele)
    return ht.group_by(ht.locus).aggregate(
        variants=hl.agg.collect(
            hl.struct(
                effect_allele=ht.effect_allele,
                noneffect_allele=ht.noneffect_allele,
                weight=ht.weight,
            )
        )
    )


def _match_weight_at_locus(
    variants_at_locus: hl.expr.ArrayExpression,
    canonical_alleles: hl.expr.ArrayExpression,
) -> hl.expr.StructExpression:
    """
    Picks the weights row at a locus whose alleles match this variant.

    `variants_at_locus` is the `variants` array from `_group_weights_by_locus`
    (missing if the locus carries no weight); `canonical_alleles` is the row's
    minimal `[ref, alt]` from `_split_multi_with_total_dosage`. Returns the
    matching `struct(effect_allele, noneffect_allele, weight)`, or missing.

    The match is on the unordered allele **set**, so it succeeds whichever way
    round the GWAS wrote the pair -- the orientation itself is resolved later,
    per row, from `canonical_alleles[0]` (the REF). `.find` on a missing array
    is missing, so a locus with no weights simply yields no match.
    """
    target = hl.set(canonical_alleles)
    return variants_at_locus.find(
        lambda v: hl.set([v.effect_allele, v.noneffect_allele]) == target
    )


def _split_multi_with_total_dosage(
    vds: hl.vds.VariantDataset,
) -> hl.MatrixTable:
    """
    Splits multi-allelic sites, preserving each entry's *total* non-ref count.

    `hl.vds.split_multi` emits one bi-allelic row per ALT allele and
    **downcodes** the genotype: at the `[C, T]` row of a `C / G,T` site, a
    sample carrying the G has its G rewritten to the reference, so its `GT` is
    `0/0` and `GT.n_alt_alleles()` is 0 -- identical to a genuine
    homozygous-reference sample.

    That is correct when counting copies of a specific ALT (the sample really
    does carry zero copies of T), but it destroys the information needed to
    count copies of the **REF**: `2 - n_alt` is only the reference count when
    `n_alt` covers *every* non-reference allele, and after downcoding it does
    not.

    So before splitting, each entry is annotated with `n_non_ref` -- the number
    of non-reference alleles in its **pre-split** local genotype, across all
    ALTs at the site. Extra entry fields survive `split_multi`, so this is
    carried through to every row the site is split into, and
    `_entry_contribution` uses it for REF-effect variants. See
    `_orient_weight_and_offset` for the arithmetic.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        An interval-filtered VariantDataset, with local (`LGT`) or global (`GT`)
        genotypes.

    Returns
    -------
    hail.MatrixTable
        The split `variant_data`, still keyed on its own (possibly non-minimal)
        `(locus, alleles)`, with a `GT` entry field (downcoded, per-ALT), an
        `n_non_ref` entry field (total, pre-split), and a `canonical_alleles`
        row field holding the minimal `[ref, alt]` the join matches on.
    """
    vd = vds.variant_data

    # All of Us ships local genotypes (LGT/LA). `LGT.n_alt_alleles()` counts
    # local non-zero indices, and local index 0 is always the reference, so
    # this is the total number of non-reference alleles the sample carries at
    # the site -- exactly what downcoding is about to destroy.
    genotype = vd.LGT if "LGT" in vd.entry.dtype else vd.GT
    vd = vd.annotate_entries(n_non_ref=genotype.n_alt_alleles())

    # `filter_changed_loci` is left at its default of False, which *raises* if
    # `hl.min_rep` moves a variant's locus. That is deliberate; do not set it
    # to True to silence an error.
    #
    # min_rep trims shared bases. Trimming a shared SUFFIX is safe, and is
    # relied upon: it reduces [AGGGC, A, GGGGC] to A/G at the same locus, which
    # is how a GWAS names that SNP. Trimming a shared PREFIX instead MOVES the
    # locus ([GG, G, GT] -> G/T one base downstream), and hail will not
    # relocate a row -- it can only raise, or drop the allele. Dropping it
    # would mean the variant silently never scores.
    #
    # No variant of the prefix-trimming shape exists in All of Us: 0 of
    # 6,001,424 ALT alleles in a 10Mb window, 21% of whose rows were
    # multi-allelic (notebooks/measure_minrep_locus_shift.ipynb). So this is
    # not a crash risk, it is a tripwire: if a future VDS release changes
    # variant representation we fail loudly, instead of scoring silently
    # wrong. Pinned by
    # tests/integration/test_allele_matching.py
    # ::test_a_locus_shifting_variant_raises_rather_than_vanishing
    split = hl.vds.split_multi(hl.vds.VariantDataset(vds.reference_data, vd))
    svd = split.variant_data

    # split_multi only guarantees minimal representation for rows it actually
    # SPLIT. An already-biallelic row is passed through with its ORIGINAL
    # alleles, so a non-minimal biallelic variant -- e.g. AAAG/GAAG, which a
    # GWAS names A/G -- keeps AAAG/GAAG and would never join the minimally-named
    # weights: it would score 0 with no error (`n_matched` just runs low). The
    # min_rep inside split_multi, and its filter_changed_loci tripwire, reach
    # only the rows it split.
    #
    # So annotate every row with its minimal alleles here. The join matches on
    # this field, not on the row key, so we do NOT re-key -- re-keying would
    # shuffle the whole chunk (entries and all); annotating is a free row op.
    # min_rep is idempotent on the rows split_multi already reduced, so this
    # only changes the biallelic passthroughs. See `_group_weights_by_locus`
    # and `_match_weight_at_locus` for the shuffle-free join that consumes it.
    #
    # The same locus-shift tripwire applies: a shared PREFIX moves the locus,
    # which hail cannot key across, so raise rather than silently drop -- the
    # biallelic analogue of split_multi's filter_changed_loci=False, which does
    # not cover these passthrough rows.
    mr = hl.min_rep(svd.locus, svd.alleles)
    canonical_alleles = (
        hl.case()
        .when(mr.locus == svd.locus, mr.alleles)
        .or_error(
            "min_rep moved a variant's locus at "
            + hl.str(svd.locus)
            + " -- a shared-prefix representation hail cannot key across. "
            "See _split_multi_with_total_dosage; do not silence this."
        )
    )
    return svd.annotate_rows(canonical_alleles=canonical_alleles)


def _orient_weight_and_offset(
    ref_allele: hl.expr.StringExpression,
    weights_info: hl.expr.StructExpression,
) -> tuple[
    hl.expr.Float64Expression,
    hl.expr.Float64Expression,
    hl.expr.BooleanExpression,
]:
    """
    Resolves a matched weights row against the VDS's own REF/ALT orientation.

    Returns `(weight_per_alt_copy, hom_ref_offset, ref_is_effect)`. The score
    for a variant is

        agg.sum(weight_per_alt_copy * dosage)  +  hom_ref_offset

    where the first term is aggregated over *entries* and the second is a
    row-level constant added to every sample. `ref_is_effect` selects which
    dosage the entry term counts; see `_entry_contribution`.

    Why the offset exists
    ---------------------
    When the effect allele is the ALT, a sample's contribution is `w * n_alt`
    and a homozygous-reference sample correctly contributes 0.

    When the effect allele is the **REF**, the contribution is `w * (2 -
    n_non_ref)`, and a hom-ref sample should get the full `2w`. But a hom-ref
    sample has **no entry** in `variant_data` -- hail filters absent entries out
    of the entry stream, so no aggregator ever visits it and no per-entry
    default can reach it. Confirmed on the real All of Us VDS: 94 of 200 samples
    had zero entries across a 5-variant window.

    Rewriting the contribution as

        w * (2 - n_non_ref)  ==  2w  -  w * n_non_ref

    splits it into a part that is already correct for absent samples (their
    `n_non_ref` is 0, so `-w * n_non_ref` is 0) and a constant `2w` that does
    not depend on the genotype at all. The constant is therefore recoverable as
    a row-level term, with no densification and no extra pass over the VDS.

    **`n_non_ref` is not `GT.n_alt_alleles()`** at a multi-allelic site. See
    `_split_multi_with_total_dosage`: after splitting, `GT` counts copies of one
    specific ALT, and a sample carrying a *different* ALT downcodes to `0/0`.
    Using it here would score that sample as two copies of the reference when it
    carries one, or none.

    The caller must only apply the offset for rows that actually matched a
    variant in the VDS. Crediting `2w` for a variant that is not in the callset
    would invent signal for every sample.

    Parameters
    ----------
    ref_allele : hail.expr.StringExpression
        The variant's reference allele in minimal form -- `canonical_alleles[0]`
        from `_split_multi_with_total_dosage`, not the row's own (possibly
        non-minimal) `alleles[0]`, so orientation lines up with the join.
    weights_info : hail.expr.StructExpression
        The matched weights row: 'effect_allele', 'noneffect_allele', 'weight'.

    Returns
    -------
    tuple
        The per-copy weight, the hom-ref offset, and whether the effect allele
        is the reference base for this row.
    """
    ref_is_effect = weights_info.effect_allele == ref_allele
    weight = weights_info.weight

    weight_per_alt_copy = hl.if_else(ref_is_effect, -weight, weight)
    hom_ref_offset = hl.if_else(ref_is_effect, 2.0 * weight, 0.0)

    return weight_per_alt_copy, hom_ref_offset, ref_is_effect


def _entry_contribution(
    mt: hl.MatrixTable,
    weight_per_alt_copy: hl.expr.Float64Expression,
    hom_ref_offset: hl.expr.Float64Expression,
    ref_is_effect: hl.expr.BooleanExpression,
) -> hl.expr.Float64Expression:
    """
    The per-entry term of the score, paired with a row-level `hom_ref_offset`.

    Which dosage is counted depends on the orientation of the row, and the two
    are **not interchangeable at a multi-allelic site**:

    * effect allele is the **ALT** -- count `GT.n_alt_alleles()`, copies of this
      row's specific ALT. Downcoding is exactly right here: a sample carrying a
      different ALT carries zero copies of this one.
    * effect allele is the **REF** -- count `n_non_ref`, the sample's *total*
      number of non-reference alleles at the site, taken from its pre-split
      genotype by `_split_multi_with_total_dosage`. Copies of the reference are
      `2 - n_non_ref`, and that identity only holds if every non-reference
      allele is counted. Using the downcoded `GT.n_alt_alleles()` would score a
      sample carrying a *different* ALT as homozygous reference -- crediting it
      two copies of an allele it holds one of, or none of.

    For a **no-call** -- an entry that exists but whose genotype is unknown --
    the contribution is `-hom_ref_offset`, which exactly cancels the offset that
    row adds to every sample. The sample therefore contributes nothing for that
    variant.

    That cancellation is the whole reason this function exists. The offset is
    added unconditionally to reach homozygous-reference samples, who have **no
    entry** and cannot be reached any other way. But a no-call *does* have an
    entry, so it receives the offset too -- and without this correction an
    unknown genotype would be scored as two copies of the reference. That is
    exactly the bug that sank the old non-split path, reappearing by a different
    route.

    Cancelling also keeps missingness handled consistently. A no-call at an
    ALT-effect variant already contributes 0 (its dosage is missing and
    `hl.agg.sum` skips it). Without the cancellation, a no-call at a REF-effect
    variant would instead contribute `2w` -- so whether an unknown genotype was
    dropped or imputed would depend on which allele the GWAS happened to label
    the effect allele, which is arbitrary.

    Parameters
    ----------
    mt : hail.MatrixTable
        A split MatrixTable with `GT` and `n_non_ref` entry fields, as produced
        by `_split_multi_with_total_dosage`.
    weight_per_alt_copy : hail.expr.Float64Expression
        From `_orient_weight_and_offset`.
    hom_ref_offset : hail.expr.Float64Expression
        From `_orient_weight_and_offset`.
    ref_is_effect : hail.expr.BooleanExpression
        From `_orient_weight_and_offset`.

    Returns
    -------
    hail.expr.Float64Expression
        The entry's contribution to the score. Never missing.
    """
    dosage = hl.if_else(ref_is_effect, mt.n_non_ref, mt.GT.n_alt_alleles())
    return hl.if_else(
        hl.is_defined(mt.GT),
        weight_per_alt_copy * dosage,
        -hom_ref_offset,
    )


def _prepare_weights_for_chunking(
    weights_table: hl.Table, config: PRSConfig, validate_table: bool = True
) -> tuple[hl.Table, int]:
    """
    Prepares and annotates a weights table for chunked processing.

    This helper function takes a raw weights table, optionally validates it,
    and assigns each row a `chunk_id` based on the configured chunk size.
    It enables iterative PRS calculation by partitioning the input into
    manageable chunks.

    Parameters
    ----------
    weights_table : hail.Table
        Raw input weights table from the user.
    weight_col_name : str
        Column name containing effect weights.
    config : PRSConfig
        Configuration object with `chunk_size`, `log_transform_weight`, and
        `detailed_timings` settings.
    validate_table : bool, default=True
        If True, validates and preprocesses the weights table using
        `_validate_and_prepare_weights_table`. If False, assumes the table is
        already validated. This is useful for batch processing where validation
        occurs upstream.

    Returns
    -------
    tuple[hail.Table, int]
        Tuple containing:
        - A validated and annotated weights table with a `chunk_id` column.
        - The number of chunks the table was divided into.

    Raises
    ------
    ValueError
        If the input table is empty after validation.
    """
    with _log_timing(
        "Preparing and analyzing weights table", config.detailed_timings
    ):
        if validate_table:
            weights_table = _validate_and_prepare_weights_table(
                weights_table=weights_table, config=config
            )

        total_variants = weights_table.count()
        if total_variants == 0:
            raise ValueError("Weights table is empty after validation.")

        effective_chunk_size = config.chunk_size or total_variants
        n_chunks = ceil(total_variants / effective_chunk_size)
        logger.info(
            "Total variants: %d, Number of chunks: %d",
            total_variants,
            n_chunks,
        )

        # Don't use chain (hl.Table.add_index().annotate()) as it is not find
        # the idx at the annotation step due to lazy eval.
        weights_table = weights_table.add_index()
        weights_table = weights_table.annotate(
            chunk_id=hl.int(weights_table.idx / effective_chunk_size)
        )
        # Note: add_index() preserves existing keys (e.g., locus),
        # so no need to re-key explicitly

        return weights_table, n_chunks


def _create_1bp_intervals(table_chunk: hl.Table) -> hl.Table:
    """
    Creates a table of 1-base-pair (1-bp) intervals from a Hail Table keyed by
    locus.

    Useful for interval-based joins or filtering operations, where each variant
    is represented as a genomic interval spanning exactly one position
    (i.e., [locus, locus], inclusive).

    Parameters
    ----------
    table_chunk : hail.Table
        A Hail Table containing a 'locus' field of type `locus<Locus>`.

    Returns
    -------
    hail.Table
        A new Table keyed by 1-bp interval around each locus.
    """
    return table_chunk.select(
        interval=hl.interval(
            table_chunk.locus, table_chunk.locus, includes_end=True
        )
    ).key_by("interval")
