"""Allele matching and dosage, against a real hail VDS.

Every weight in the mock GWAS summary is 1.0, so a sample's `prs` is exactly the
number of effect-allele copies the library counted for it. Each expected value
below is derived by hand from the genotypes in `conftest.VARIANTS` and is
labelled either CORRECT or BUG. Tests marked BUG pin down what the code does
today so that a fix is reviewable -- they do not endorse it. See TODO.md.

The mock VDS (weights effect allele in brackets):

  locus       VDS alleles        S1   S2   S3   S4       weights   note
  chr1:1000   A / G              A/A  A/G  G/G  A/A      [G] / A
  chr1:2000   A / G              A/A  A/G  G/G  no-call  [A] / G
  chr1:3000   A / T              A/A  A/T  T/T  A/A      [A] / G   no A/G here
  chr1:4000   A / T              A/A  A/T  T/T  A/A      [G] / A   no A/G here
  chr1:5000   C / G,T            C/C  C/G  C/T  C/C      [T] / C
  chr1:6000   AGGGC / A,GGGGC    -/-  -/G  G/G  -/-      [G] / A   min_rep'd
  chr1:7000   A / C,G            A/A  C/C  A/G  A/A      [A] / G
  chr1:8000   G / A              G/G  G/A  A/A  G/G      both      REF > ALT
  chr1:8500   AAAG / GAAG        -/-  -/G  G/G  -/-      [G] / A   biallelic mr

S1 is homozygous reference at every site, so it has no entry in `variant_data`
at all. S4 has an entry at chr1:2000 whose genotype is missing -- a no-call.

chr1:8000 is the only site whose REF does not sort before its ALT. That is not
cosmetic: the weights join is keyed on alleles, and every other row here would
join just as well with the two sides keyed inconsistently. See
`test_a_variant_whose_ref_sorts_after_its_alt_scores`.

There is one scoring path: multi-allelics are split, the weights are joined on
(locus, allele pair) -- emitted in both orders, so the VDS's own REF/ALT order
is what the key matches -- and each matched row is oriented against the VDS's
own REF/ALT. A REF-effect row contributes `-w * n_non_ref` per entry plus a
row-level `2w` offset applied to every sample -- the only way to reach a hom-ref
sample, which has no entry at all.

`n_non_ref` is the sample's **total** non-reference allele count, read off its
pre-split genotype. It is deliberately *not* the post-split
`GT.n_alt_alleles()`, which counts copies of one specific ALT: splitting
downcodes every other ALT to the reference, so at a multi-allelic site a carrier
of a different ALT would otherwise be scored as homozygous reference. chr1:5000
and chr1:7000 exist to pin that.

Two things were removed. The non-split path silently dropped every hom-ref
sample at a REF-effect variant, reordering the cohort. The global
`ref_is_effect_allele` flag silently dropped every weights row of the opposite
orientation. Both are now a hard `TypeError`; see `tests/prs/test_config.py`.
"""

import hail as hl
import pytest

from aoutools.prs._calculator import _calculate_prs_chunk
from aoutools.prs._calculator_batch import (
    _calculate_prs_chunk_batch,
    _prepare_batch_weights_data,
)
from aoutools.prs._calculator_utils import _validate_and_prepare_weights_table
from aoutools.prs._config import PRSConfig

pytestmark = pytest.mark.integration


def score(vds, raw_weights, positions=None, **config_kwargs):
    """Scores the mock VDS and returns {sample: prs} plus n_matched.

    `positions` optionally restricts the weights to a subset of loci, so a
    single scenario can be isolated from the other four.
    """
    config = PRSConfig(include_n_matched=True, **config_kwargs)
    weights = _validate_and_prepare_weights_table(raw_weights, config)
    if positions is not None:
        keep = hl.set(positions)
        weights = weights.filter(keep.contains(weights.locus.position))

    df = _calculate_prs_chunk(weights, vds, config).to_pandas()
    prs = dict(zip(df["person_id"], df["prs"], strict=True))
    n_matched = int(df["n_matched"].iloc[0])
    return prs, n_matched


# --------------------------------------------------------------------------
# The fact that drives everything else.
# --------------------------------------------------------------------------


def test_hom_ref_samples_are_absent_from_the_entry_stream(vds_lgt):
    """A hom-ref sample is not a missing entry -- it is not an entry.

    Hail *filters* absent entries out of the entry stream, so the aggregators
    never visit S1 at all and no default of any kind can be applied to it. Any
    scheme that tries to give a hom-ref sample a dosage by handling a missing
    genotype is therefore unreachable -- which is what sank the non-split path,
    and what Task 2's constant-offset fix has to work around.

    Confirmed on the real All of Us VDS: 94 of 200 samples had zero entries
    across a 5-variant window. See `notebooks/verify_hom_ref_dosage.ipynb`.
    """
    vd = vds_lgt.variant_data

    visited = vd.select_cols(n=hl.agg.count()).cols().to_pandas()
    per_sample = dict(zip(visited["s"], visited["n"], strict=True))

    # Taken from the VDS, not hardcoded: adding a fixture variant must not
    # silently turn this into an assertion about the wrong number of sites.
    n_sites = vd.count_rows()

    assert per_sample["S1"] == 0, "hom-ref sample must not be visited at all"
    assert per_sample["S2"] == n_sites  # a call at every site
    assert per_sample["S3"] == n_sites
    assert per_sample["S4"] == 1  # the no-call entry at chr1:2000 exists


# --------------------------------------------------------------------------
# Dosage.
# --------------------------------------------------------------------------


def test_counts_alt_copies(vds_lgt, raw_weights):
    """CORRECT. Where the effect allele is the ALT, dosage is a plain count of
    ALT copies, the offset is zero, and a hom-ref sample contributes nothing --
    which is right, because it carries no copies of the effect allele.

    chr1:1000 (effect G) and chr1:5000 (effect T) are the two ALT-effect rows.
    """
    prs, n_matched = score(vds_lgt, raw_weights, positions=[1000, 5000])

    assert n_matched == 2
    assert prs == {
        "S1": 0.0,  # hom-ref at both: no entry, and no offset to receive
        "S2": 1.0,  # one G at 1000; C/G at 5000 carries no T
        "S3": 3.0,  # two G at 1000, one T at 5000
        "S4": 0.0,  # hom-ref at both
    }


def test_handles_multiallelic_sites(vds_lgt):
    """CORRECT. Splitting chr1:5000 (C / G,T) yields a [C,T] row that the
    weights join to, and only the T carrier scores. S2 carries the *other* ALT
    and downcodes to 0/0."""
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 5000, "effect_allele": "T",
          "noneffect_allele": "C", "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip
    prs, n_matched = score(vds_lgt, raw)

    assert n_matched == 1
    assert prs["S2"] == 0.0  # C/G -- carries the other ALT
    assert prs["S3"] == 1.0  # C/T


def test_normalizes_a_non_minimal_representation(vds_lgt):
    """CORRECT, and this is the normalization the library actually relies on.

    `chr1:6000` is stored as `[AGGGC, A, GGGGC]` -- a deletion sharing a record
    with a SNP. No GWAS names that SNP that way; it names it `A/G`. `min_rep`,
    applied by `hl.vds.split_multi`, trims the shared **suffix** GGGC and
    reduces the pair to `A/G` at this same locus -- exactly the key the weights
    join on. A multi-allelic, non-minimally-represented variant is therefore
    scored correctly.

    The contrast that matters: trimming a shared **prefix** would instead *move*
    the locus, and `split_multi` refuses to relocate rows. See
    `test_a_locus_shifting_variant_raises_rather_than_vanishing`. Suffix
    trimming is safe; prefix trimming is not. Only the former occurs in AoU.
    """
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 6000, "effect_allele": "G",
          "noneffect_allele": "A", "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip
    prs, n_matched = score(vds_lgt, raw)

    assert n_matched == 1, "the A/G weights row must match [AGGGC, A, GGGGC]"
    assert prs["S1"] == 0.0  # hom-ref
    assert prs["S2"] == 1.0  # carries one copy of the SNP allele
    assert prs["S3"] == 2.0  # homozygous for it
    assert prs["S4"] == 0.0  # hom-ref


def test_normalizes_a_non_minimal_biallelic_variant(vds_lgt):
    """CORRECT, and the case the multi-allelic test above does NOT cover.

    `chr1:8500` is stored as `[AAAG, GAAG]` -- an A->G SNP written with three
    shared suffix bases, biallelic. A GWAS names it `A/G`. Crucially,
    `hl.vds.split_multi` min_reps only the rows it actually *splits*; a row that
    is already biallelic is passed through with its original alleles. So the
    normalization the multi-allelic case gets for free (chr1:6000) does NOT
    happen here unless the library does it explicitly.

    Without that step this row stays `AAAG/GAAG`, the `A/G` weights key never
    matches, and every carrier scores 0 with no error -- exactly the FAIL that
    validate_scoring_on_aou.ipynb turned up on the real VDS (chr1:1409159).
    """
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 8500, "effect_allele": "G",
          "noneffect_allele": "A", "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip
    prs, n_matched = score(vds_lgt, raw)

    assert n_matched == 1, (
        "the A/G weights row must match the biallelic [AAAG, GAAG]; if it does "
        "not, split_multi left the passthrough row non-minimal"
    )
    assert prs["S1"] == 0.0  # hom-ref
    assert prs["S2"] == 1.0  # carries one copy of the SNP allele
    assert prs["S3"] == 2.0  # homozygous for it
    assert prs["S4"] == 0.0  # hom-ref


def test_ref_effect_at_a_multiallelic_site(vds_lgt):
    """CORRECT, and the reason `n_non_ref` exists. This was Finding 6.

    `hl.vds.split_multi` **downcodes**: at the `[C, T]` row of the `C / G,T`
    site, a sample carrying the G has its G rewritten to the reference. Its `GT`
    becomes `0/0` and `GT.n_alt_alleles()` is 0 -- byte-for-byte identical to a
    genuine homozygous-reference sample.

    For an ALT-effect weight that is exactly right: S2 really does carry zero
    copies of T. But for a **REF-effect** weight it is a disaster. The row-level
    offset hands every sample `2w`, and the entry term subtracts `w * dosage`
    back off. If `dosage` is the downcoded count, S2 subtracts nothing and is
    scored as carrying **two** copies of C -- when it carries one.

    The fix is `n_non_ref`: the sample's total non-reference allele count, read
    off its **pre-split** local genotype and carried through the split by
    `_split_multi_with_total_dosage`. Copies of the reference are
    `2 - n_non_ref`, and that identity only holds when *every* non-reference
    allele is counted.

    chr1:5000 is C / G,T. Effect allele C (the REF), paired with T.
    True copies of C: S1 2 (hom-ref), S2 1 (C/G), S3 1 (C/T), S4 2 (hom-ref).
    """
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 5000, "effect_allele": "C",
          "noneffect_allele": "T", "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip
    prs, n_matched = score(vds_lgt, raw)

    assert n_matched == 1
    assert prs["S1"] == 2.0  # C/C, no entry: the offset is the whole score
    assert prs["S2"] == 1.0, (
        "S2 is C/G. It carries ONE C. Its G downcodes to REF at the [C,T] row, "
        "so a dosage read from the downcoded GT would score it 2.0."
    )
    assert prs["S3"] == 1.0  # C/T
    assert prs["S4"] == 2.0  # C/C, no entry


def test_ref_effect_when_a_sample_is_homozygous_for_another_alt(vds_lgt):
    """CORRECT. The same bug as above, at its worst.

    chr1:7000 is A / C,G and S2 is **C/C** -- homozygous for an ALT the weights
    row never names. At the split `[A, G]` row its genotype downcodes to `0/0`,
    so a downcoded dosage would credit it the full `2w`: two copies of A, an
    allele it does not carry at all. The error is a whole `2w`, the largest a
    single variant can produce.

    Effect allele A (the REF), paired with G.
    True copies of A: S1 2 (hom-ref), S2 0 (C/C), S3 1 (A/G), S4 2 (hom-ref).
    """
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 7000, "effect_allele": "A",
          "noneffect_allele": "G", "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip
    prs, n_matched = score(vds_lgt, raw)

    assert n_matched == 1
    assert prs["S1"] == 2.0  # A/A, no entry
    assert prs["S2"] == 0.0, (
        "S2 is C/C and carries NO copy of A. A dosage read from the downcoded "
        "GT would score it 2.0 -- the maximum possible error for one variant."
    )
    assert prs["S3"] == 1.0  # A/G
    assert prs["S4"] == 2.0  # A/A, no entry


@pytest.mark.parametrize(
    ("effect", "noneffect", "expected"),
    [
        # Effect allele is the ALT ('A'): copies of A are 0, 1, 2, 0.
        ("A", "G", {"S1": 0.0, "S2": 1.0, "S3": 2.0, "S4": 0.0}),
        # Effect allele is the REF ('G'): copies of G are 2, 1, 0, 2. This one
        # also needs the hom-ref offset, so it fails differently if the join is
        # right but the orientation is read off the wrong allele.
        ("G", "A", {"S1": 2.0, "S2": 1.0, "S3": 0.0, "S4": 2.0}),
    ],
)
def test_a_variant_whose_ref_sorts_after_its_alt_scores(
    vds_lgt, effect, noneffect, expected
):
    """CORRECT, and the regression guard for a silent half-the-file dropout.

    chr1:8000 is `G / A`: its REF sorts **after** its ALT. Every other variant
    in the mock VDS has REF < ALT, and that accident hid a real bug for the
    whole life of the suite.

    The weights used to be keyed on the *sorted* allele pair (`[A, G]`), but
    the split `variant_data` row key is `(locus, [REF, ALT])` in the VDS's own
    order (`[G, A]`). The two never met. The variant simply did not join, and a
    weights row that does not join is indistinguishable from one whose variant
    is absent from the callset: no error, no warning, a finite score.

    On a real PGS this dropped ~half of every file -- 922 of 1,940 variants for
    PGS000746 against the All of Us VDS -- and biased what survived toward
    REF=A and REF=C. Both orientations are checked here because the fix
    duplicates the weights row rather than sorting the VDS key, and getting
    that backwards would repair the join while breaking `alleles[0] == REF`,
    which is what `_orient_weight_and_offset` reads.

    Genotypes are identical to chr1:1000, so the expectations are the same copy
    counts, just written against the mirrored allele order.
    """
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 8000, "effect_allele": effect,
          "noneffect_allele": noneffect, "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip
    prs, n_matched = score(vds_lgt, raw)

    assert n_matched == 1, (
        f"the G/A variant at chr1:8000 did not join at all (effect={effect}). "
        "Its REF sorts after its ALT -- if the weights key is sorted but the "
        "VDS row key is not, this variant silently vanishes from the score."
    )
    assert prs == expected


def test_alt_effect_is_unaffected_by_another_alt(vds_lgt):
    """CORRECT, and the other half of the argument.

    Downcoding is only wrong for REF-effect rows. When the effect allele is an
    ALT, a sample carrying a *different* ALT genuinely carries zero copies of
    the effect allele, and the downcoded `GT.n_alt_alleles()` says exactly that.

    So the two orientations must count *different* dosages, which is why
    `_entry_contribution` branches on `ref_is_effect` rather than using one
    dosage everywhere. Counting `n_non_ref` here instead would score S2 -- a
    C/G carrier -- as holding one copy of T.

    chr1:5000 is C / G,T. Effect allele T (an ALT).
    True copies of T: S1 0, S2 0 (C/G), S3 1 (C/T), S4 0.
    """
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 5000, "effect_allele": "T",
          "noneffect_allele": "C", "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip
    prs, _ = score(vds_lgt, raw)

    assert prs == {"S1": 0.0, "S2": 0.0, "S3": 1.0, "S4": 0.0}


def test_a_locus_shifting_variant_raises_rather_than_vanishing(
    vds_locus_shifting,
):
    """A TRIPWIRE, deliberately left armed. Do not "fix" this by passing
    `filter_changed_loci=True`.

    `chr1:1001 [GG, G, GT]` minreps its SNP allele to `G/T` at chr1:**1002**.
    `hl.vds.split_multi` will not relocate a row, so it offers only two
    behaviours: raise (`filter_changed_loci=False`, the default) or silently
    drop the allele (`True`). Neither one scores the variant.

    We keep the raising default on purpose. **No variant of this shape exists
    in All of Us**: 0 of 6,001,424 ALT alleles in a 10Mb window of chr1
    shifted, and 21% of those rows were multi-allelic
    (`notebooks/measure_minrep_locus_shift.ipynb`). So the exception is not a
    crash risk -- it is a tripwire. If a future VDS release ever changes the
    variant representation, the run fails loudly instead of quietly dropping
    variants and producing a plausible, wrong score.

    Setting `filter_changed_loci=True` would convert that tripwire into exactly
    the kind of silent data loss the rest of this file exists to prevent.
    """
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 1002, "effect_allele": "T",
          "noneffect_allele": "G", "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip

    with pytest.raises(Exception, match="non-left-aligned"):
        score(vds_locus_shifting, raw)


def test_does_not_invent_a_genotype_for_a_no_call(vds_lgt, raw_weights):
    """CORRECT, and it is the subtlest guard in the file.

    S4's entry at chr1:2000 exists but its genotype is missing -- a genuine
    no-call, an *unknown* genotype. The effect allele there is A, the REF, so
    that row carries a `hom_ref_offset` of 2w which is added to **every** sample
    in order to reach the hom-ref samples who have no entry.

    S4 would therefore collect the 2w as well, and an unknown genotype would be
    scored as two copies of the reference -- the exact bug that sank the old
    non-split path, arriving by a new route. `_entry_contribution` cancels it:
    a no-call has an entry, so it can be given `-hom_ref_offset`, which zeroes
    out. S4 contributes nothing for this variant, as it should.

    Compare S1, who is hom-ref here and has *no* entry to cancel with, and so
    correctly keeps the full 2w.
    """
    prs, _ = score(vds_lgt, raw_weights, positions=[2000])

    assert prs["S4"] == 0.0, "an unknown genotype must not be scored as A/A"
    assert prs["S1"] == 2.0, "a hom-ref sample must still receive the offset"


# --------------------------------------------------------------------------
# Allele orientation, resolved per row against the VDS.
# --------------------------------------------------------------------------


def test_scores_a_ref_effect_variant_exactly(vds_lgt, raw_weights):
    """CORRECT, and this is the test the whole hom-ref problem comes down to.

    At chr1:2000 the effect allele is A -- the REF base. True copies of A:
    S1 A/A = 2, S2 A/G = 1, S3 G/G = 0. S1 has no entry in `variant_data` and is
    never visited by any entry aggregator, so its 2 copies cannot come from the
    dosage. They arrive as the row-level `hom_ref_offset` (2w), which is added
    to every sample.

    These are the truth, not an offset from it. Before this fix the scores were
    0/-1/-2 -- correctly *ordered*, but shifted by a constant 2w, so absolute
    values were meaningless.
    """
    prs, n_matched = score(vds_lgt, raw_weights, positions=[2000])

    assert n_matched == 1
    assert prs["S1"] == 2.0  # A/A -- supplied entirely by the offset
    assert prs["S2"] == 1.0  # A/G
    assert prs["S3"] == 0.0  # G/G
    # S4 is a no-call; its offset is cancelled. See
    # test_does_not_invent_a_genotype_for_a_no_call.
    assert prs["S4"] == 0.0


def test_scores_both_orientations_from_one_file(vds_lgt, raw_weights):
    """CORRECT. Orientation is a per-row property, read off the VDS.

    The mock weights file is mixed: chr1:1000 and chr1:5000 name the ALT as
    the effect allele, chr1:2000 names the REF. The old global
    `ref_is_effect_allele` flag could score one group or the other -- 2 rows
    with it off, 1 with it on -- and silently dropped the rest, not even
    counting them in `n_matched`. The PGS Catalog does not harmonize the effect
    allele onto the ALT, so a mixed file is the normal case.

    All three real variants now match in a single pass.
    """
    prs, n_matched = score(vds_lgt, raw_weights)

    assert n_matched == 3, "1000 and 5000 (ALT-effect) plus 2000 (REF-effect)"
    assert prs == {
        # 1000: one G. 2000: one A. 5000: C/G, no T.
        "S2": 1.0 + 1.0 + 0.0,
        # 1000: two G. 2000: no A. 5000: one T.
        "S3": 2.0 + 0.0 + 1.0,
        # hom-ref everywhere: 0 at 1000 and 5000, but two copies of A at 2000.
        "S1": 0.0 + 2.0 + 0.0,
        # hom-ref at 1000 and 5000; its 2000 entry is a no-call, contributing
        # nothing rather than being imputed to A/A.
        "S4": 0.0 + 0.0 + 0.0,
    }


def test_does_not_credit_the_offset_for_an_unmatched_variant(
    vds_lgt, raw_weights
):
    """CORRECT. Variant identity, and the guard on the offset -- one test,
    because they share an assertion.

    chr1:3000 names an A/G SNP with A (a REF base) as the effect allele, and
    chr1:4000 an A/G SNP with G as the effect allele. The VDS has an **A/T**
    SNP at both positions; the A/G variant does not exist. Keying the join on
    (locus, allele pair) means neither row can match -- that key *is* the
    variant-identity check, and it is what the removed path had to opt into via
    `strict_allele_match`.

    Crucially the 3000 row must also contribute **no offset**. Its effect
    allele is a REF base, so a matched row would carry a `2w` credit for every
    sample. Handing that out for a variant that is not in the callset would
    invent signal for the entire cohort, including samples with no entry
    anywhere.

    The single-score path gets this from `filter_rows` (the row is gone before
    the offset is aggregated). The batch path never filters rows, so it must
    gate the offset on `is_valid_{score}` -- a different mechanism for the same
    guarantee, which is why `test_batch_agrees_with_single` matters.
    """
    prs, n_matched = score(vds_lgt, raw_weights, positions=[3000, 4000])

    assert n_matched == 0, "neither A/G variant exists in the VDS"
    assert set(prs.values()) == {0.0}, (
        "an unmatched weights row must contribute nothing at all -- not even "
        f"the hom-ref offset; got {prs}"
    )


# --------------------------------------------------------------------------
# Batch must agree with single.
# --------------------------------------------------------------------------


def test_batch_agrees_with_single(vds_lgt, raw_weights):
    """The batch calculator masks with `hl.if_else` where the single-score path
    filters rows. Different mechanism, same number required -- including for the
    hom-ref offset, which batch has to gate explicitly."""
    config = PRSConfig(include_n_matched=True)
    single, _ = score(vds_lgt, raw_weights)

    prepared, _ = _prepare_batch_weights_data({"s1": raw_weights}, config)
    df = _calculate_prs_chunk_batch(
        vds_lgt, {"s1": raw_weights}, prepared, config
    ).to_pandas()
    batch = dict(zip(df["s"], df["s1"], strict=True))

    assert batch == single
