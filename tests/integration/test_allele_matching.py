"""Allele matching and dosage, against a real hail VDS.

Every weight in the mock GWAS summary is 1.0, so a sample's `prs` is exactly the
number of effect-allele copies the library counted for it. Each expected value
below is derived by hand from the genotypes in `conftest.VARIANTS` and is
labelled either CORRECT or BUG. Tests marked BUG pin down what the code does
today so that a fix is reviewable -- they do not endorse it. See TODO.md.

The mock VDS (weights effect allele in brackets):

  locus       VDS alleles  S1   S2   S3   S4       weights   note
  chr1:1000   A / G        A/A  A/G  G/G  A/A      [G] / A
  chr1:2000   A / G        A/A  A/G  G/G  no-call  [A] / G
  chr1:3000   A / T        A/A  A/T  T/T  A/A      [A] / G   no A/G here
  chr1:4000   A / T        A/A  A/T  T/T  A/A      [G] / A   no A/G here
  chr1:5000   C / G,T      C/C  C/G  C/T  C/C      [T] / C

S1 is homozygous reference at every site, so it has no entry in `variant_data`
at all. S4 has an entry at chr1:2000 whose genotype is missing -- a no-call.
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

BOTH_ENCODINGS = pytest.mark.parametrize("vds_fixture", ["vds_lgt", "vds_gt"])


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
# The finding that drives everything below.
# --------------------------------------------------------------------------


@BOTH_ENCODINGS
def test_hom_ref_samples_are_absent_from_the_entry_stream(vds_fixture, request):
    """A hom-ref sample is not a missing entry -- it is not an entry.

    `_calculate_dosage` has a branch that treats a missing genotype as
    homozygous reference, documented as "accounting for sparse storage of
    homozygous reference calls". This test shows that branch cannot be reached
    that way: hail *filters* absent entries out of the entry stream, so the
    aggregators never visit S1 at all and no default of any kind is applied.

    Everything surprising in this file follows from this one fact.
    """
    vds = request.getfixturevalue(vds_fixture)
    vd = vds.variant_data

    visited = vd.select_cols(n=hl.agg.count()).cols().to_pandas()
    per_sample = dict(zip(visited["s"], visited["n"], strict=True))

    assert per_sample["S1"] == 0, "hom-ref sample must not be visited at all"
    assert per_sample["S2"] == 5  # a call at every one of the five sites
    assert per_sample["S3"] == 5
    assert per_sample["S4"] == 1  # the no-call entry at chr1:2000 exists


# --------------------------------------------------------------------------
# Split path -- the default.
# --------------------------------------------------------------------------


def test_split_default_counts_alt_copies(vds_lgt, raw_weights):
    """CORRECT. Dosage is a count of ALT copies, hom-ref contributes nothing."""
    prs, n_matched = score(
        vds_lgt, raw_weights, split_multi=True, ref_is_effect_allele=False
    )

    # Only two weights rows have their effect allele on the ALT side and a
    # variant that actually exists: chr1:1000 (G) and chr1:5000 (T).
    assert n_matched == 2
    assert prs == {
        "S1": 0.0,  # hom-ref everywhere
        "S2": 1.0,  # one G at 1000; C/G at 5000 carries no T
        "S3": 3.0,  # two G at 1000, one T at 5000
        "S4": 0.0,  # its only entry is the no-call at 2000, which is unmatched
    }


def test_split_silently_drops_rows_of_the_opposite_orientation(
    vds_lgt, raw_weights
):
    """BUG (design limit). `ref_is_effect_allele` is global, not per-row.

    The split path builds the join key from the flag, so a weights file with
    mixed orientation -- some rows naming the REF as the effect allele, some the
    ALT -- has whichever half disagrees with the flag silently dropped. Three of
    the five rows here vanish with the flag off, one with it on, and the score
    still comes out as a clean number.
    """
    _, n_alt_oriented = score(
        vds_lgt, raw_weights, split_multi=True, ref_is_effect_allele=False
    )
    _, n_ref_oriented = score(
        vds_lgt, raw_weights, split_multi=True, ref_is_effect_allele=True
    )

    assert n_alt_oriented == 2  # 1000, 5000
    assert n_ref_oriented == 1  # 2000 only
    # No setting of the flag scores all five rows.
    assert n_alt_oriented + n_ref_oriented < len(raw_weights.collect())


def test_split_ref_is_effect_offsets_every_sample_equally(vds_lgt, raw_weights):
    """CORRECT, with the documented caveat.

    With `ref_is_effect_allele=True` the code negates the weight and multiplies
    by the ALT count, giving `-w * n_alt` where the truth is `w * (2 - n_alt)`.
    The two differ by a constant `2w` per matched variant. This test pins the
    thing that actually matters: the offset is the *same for every sample*,
    including the hom-ref sample that is never visited. Absolute scores are
    shifted; rankings are not. If a change ever makes this offset vary per
    sample, cross-sample comparability is destroyed -- and this test goes red.
    """
    # chr1:2000 alone. Effect allele is A (the REF). True copies of A:
    # S1 A/A = 2, S2 A/G = 1, S3 G/G = 0.
    truth = {"S1": 2.0, "S2": 1.0, "S3": 0.0}

    prs, n_matched = score(
        vds_lgt,
        raw_weights,
        positions=[2000],
        split_multi=True,
        ref_is_effect_allele=True,
    )
    assert n_matched == 1
    assert prs["S1"] == 0.0 and prs["S2"] == -1.0 and prs["S3"] == -2.0

    offsets = {s: truth[s] - prs[s] for s in truth}
    assert set(offsets.values()) == {2.0}, (
        "the 2w offset must be identical for every sample, or scores stop "
        f"being comparable across individuals; got {offsets}"
    )


# --------------------------------------------------------------------------
# Non-split path -- allele matching.
# --------------------------------------------------------------------------


@BOTH_ENCODINGS
def test_non_split_strict_drops_variants_that_are_not_there(
    vds_fixture, request, raw_weights
):
    """CORRECT. The A/G weights rows at 3000/4000 have no A/G variant in the
    VDS, and strict matching drops them for everyone."""
    vds = request.getfixturevalue(vds_fixture)
    prs, n_matched = score(
        vds, raw_weights, split_multi=False, strict_allele_match=True
    )

    # 1000, 2000, 5000 match; 3000 and 4000 do not.
    assert n_matched == 3
    assert prs == {
        "S1": 0.0,  # BUG, see test_non_split_ref_effect_loses_hom_ref_samples
        "S2": 2.0,  # one G at 1000, one A at 2000, no T at 5000
        "S3": 3.0,  # two G at 1000, no A at 2000, one T at 5000
        "S4": 2.0,  # BUG, see test_non_split_scores_a_no_call_as_hom_ref
    }


@BOTH_ENCODINGS
def test_non_split_loose_scores_a_variant_the_gwas_never_studied(
    vds_fixture, request, raw_weights
):
    """BUG. `strict_allele_match=False` performs no allele check at all.

    At chr1:3000 the weights describe an A/G SNP. The VDS has an A/T SNP. The
    join is on locus only, so the row matches; the effect allele A is the REF
    base, which is the same string whichever ALT sits there, so copies of A are
    counted and scored. The counts are *truthful* counts of A -- the defect is
    that they are weighted by an effect size estimated for an A-vs-G contrast
    while varying with each sample's unrelated T genotype.

    Concretely: S2 alone picks up a phantom copy, purely because it happens to
    be A/T rather than T/T. That is per-sample noise, so it moves rankings.
    """
    vds = request.getfixturevalue(vds_fixture)
    strict, _ = score(
        vds, raw_weights, split_multi=False, strict_allele_match=True
    )
    loose, n_matched = score(
        vds, raw_weights, split_multi=False, strict_allele_match=False
    )

    assert n_matched == 5  # every row "matched", including the two phantoms
    assert loose["S2"] - strict["S2"] == 1.0, "phantom copy of A at chr1:3000"
    # The others are unchanged only by luck: S3 is T/T so it has no A to count,
    # and S1 is never visited.
    assert loose["S3"] == strict["S3"]
    assert loose["S1"] == strict["S1"]


@BOTH_ENCODINGS
def test_non_split_loose_fails_safe_when_the_effect_allele_is_absent(
    vds_fixture, request, raw_weights
):
    """CORRECT-ish. The mirror of the case above fails safe -- but is still
    counted as a match.

    At chr1:4000 the effect allele is G, which is not in the VDS's A/T variant.
    Nothing string-matches G, so dosage is 0 for everyone and the row cannot
    corrupt a score. It is nonetheless reported in `n_matched`, which therefore
    over-reports overlap with the VDS.
    """
    vds = request.getfixturevalue(vds_fixture)
    prs, n_matched = score(
        vds,
        raw_weights,
        positions=[4000],
        split_multi=False,
        strict_allele_match=False,
    )

    assert set(prs.values()) == {0.0}, "no sample carries G at chr1:4000"
    assert n_matched == 1, "counted as matched despite contributing nothing"


@BOTH_ENCODINGS
def test_non_split_handles_multiallelic_sites(vds_fixture, request):
    """CORRECT. Local (LA) and global (GT) encodings both decode to C/G and
    C/T, and only the T carrier scores."""
    vds = request.getfixturevalue(vds_fixture)
    raw = hl.Table.parallelize(
        [{"chr": "chr1", "pos": 5000, "effect_allele": "T",
          "noneffect_allele": "C", "weight": 1.0}],
        hl.tstruct(
            chr=hl.tstr, pos=hl.tint32, effect_allele=hl.tstr,
            noneffect_allele=hl.tstr, weight=hl.tfloat64,
        ),
    )  # fmt: skip
    prs, n_matched = score(
        vds, raw, split_multi=False, strict_allele_match=True
    )

    assert n_matched == 1
    assert prs["S2"] == 0.0  # C/G -- carries the other ALT
    assert prs["S3"] == 1.0  # C/T


# --------------------------------------------------------------------------
# The two ways the non-split path mishandles a missing genotype.
# --------------------------------------------------------------------------


@BOTH_ENCODINGS
def test_non_split_ref_effect_loses_hom_ref_samples(
    vds_fixture, request, raw_weights
):
    """BUG, and the severe one: it changes rankings, not just absolute scores.

    At chr1:2000 the effect allele is A, the REF. The true copy counts are
    S1 A/A = 2, S2 A/G = 1, S3 G/G = 0. `_calculate_dosage` intends to supply
    the 2 for S1 via its missing-genotype branch, but S1 has no entry, so that
    branch is never evaluated and S1 contributes nothing.

    S2 and S3 are scored correctly. S1 is not. The error is therefore *not* a
    uniform offset like the split path's -- it applies only to samples who are
    hom-ref at the site, which is genotype-dependent and per-sample. The true
    ordering S1 > S2 > S3 comes out as S2 > S1 = S3.
    """
    vds = request.getfixturevalue(vds_fixture)
    prs, _ = score(
        vds,
        raw_weights,
        positions=[2000],
        split_multi=False,
        strict_allele_match=True,
    )

    assert prs["S2"] == 1.0  # correct
    assert prs["S3"] == 0.0  # correct
    assert prs["S1"] == 0.0, "should be 2.0 -- S1 is A/A and A is the effect"

    # The ranking inversion, stated outright.
    assert prs["S1"] < prs["S2"], (
        "S1 truly carries more A than S2, yet scores lower"
    )
    assert prs["S1"] == prs["S3"], "S1 (two copies) is tied with S3 (zero)"


@BOTH_ENCODINGS
def test_non_split_scores_a_no_call_as_hom_ref(
    vds_fixture, request, raw_weights
):
    """BUG. The missing-genotype branch fires for exactly the wrong sample.

    S4's entry at chr1:2000 exists but its genotype is missing -- a genuine
    no-call, an *unknown* genotype. Because the entry exists, it is visited, the
    missing-genotype branch fires, and S4 is handed two copies of the effect
    allele it may not have.

    So the branch that was meant to serve hom-ref samples never sees them, and
    instead invents data for no-calls. It is exactly inverted.
    """
    vds = request.getfixturevalue(vds_fixture)
    prs, _ = score(
        vds,
        raw_weights,
        positions=[2000],
        split_multi=False,
        strict_allele_match=True,
    )

    assert prs["S4"] == 2.0, "an unknown genotype is scored as two copies of A"


# --------------------------------------------------------------------------
# Batch must agree with single.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("strict", [True, False])
def test_batch_agrees_with_single(vds_lgt, raw_weights, strict):
    """The batch calculator gates matches with `hl.if_else` rather than
    `filter_rows`. Different mechanism, same number required."""
    config = PRSConfig(
        split_multi=False, strict_allele_match=strict, include_n_matched=True
    )
    single, _ = score(
        vds_lgt, raw_weights, split_multi=False, strict_allele_match=strict
    )

    prepared, _ = _prepare_batch_weights_data({"s1": raw_weights}, config)
    df = _calculate_prs_chunk_batch(
        vds_lgt, {"s1": raw_weights}, prepared, config
    ).to_pandas()
    batch = dict(zip(df["s"], df["s1"], strict=True))

    assert batch == single
