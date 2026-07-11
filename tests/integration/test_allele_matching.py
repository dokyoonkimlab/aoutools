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

There is one scoring path: multi-allelics are split and the weights are joined
on (locus, alleles). The non-split path was removed -- it silently dropped every
hom-ref sample at a REF-effect variant, which reordered the cohort. The tests
that pinned that bug are gone with it; `tests/prs/test_config.py` asserts the
config knobs are now a hard error.
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

    assert per_sample["S1"] == 0, "hom-ref sample must not be visited at all"
    assert per_sample["S2"] == 5  # a call at every one of the five sites
    assert per_sample["S3"] == 5
    assert per_sample["S4"] == 1  # the no-call entry at chr1:2000 exists


# --------------------------------------------------------------------------
# Dosage.
# --------------------------------------------------------------------------


def test_counts_alt_copies(vds_lgt, raw_weights):
    """CORRECT. Dosage is a count of ALT copies, hom-ref contributes nothing."""
    prs, n_matched = score(vds_lgt, raw_weights)

    # Only two weights rows have their effect allele on the ALT side and a
    # variant that actually exists: chr1:1000 (G) and chr1:5000 (T).
    assert n_matched == 2
    assert prs == {
        "S1": 0.0,  # hom-ref everywhere
        "S2": 1.0,  # one G at 1000; C/G at 5000 carries no T
        "S3": 3.0,  # two G at 1000, one T at 5000
        "S4": 0.0,  # its only entry is the no-call at 2000, which is unmatched
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


def test_does_not_invent_a_genotype_for_a_no_call(vds_lgt, raw_weights):
    """CORRECT, and it is a regression guard.

    S4's entry at chr1:2000 exists but its genotype is missing -- a genuine
    no-call, an *unknown* genotype. The deleted `_calculate_dosage` handed such
    a sample two copies of the effect allele whenever that allele was the REF
    base. `GT.n_alt_alleles()` is missing instead, and `hl.agg.sum` skips it, so
    S4 simply does not contribute. Nothing is invented.
    """
    prs, _ = score(
        vds_lgt, raw_weights, positions=[2000], ref_is_effect_allele=True
    )

    assert prs["S4"] == 0.0, "an unknown genotype must not be scored as A/A"


# --------------------------------------------------------------------------
# Variant identity: the join key is the allele check.
# --------------------------------------------------------------------------


def test_drops_variants_that_are_not_there(vds_lgt, raw_weights):
    """CORRECT. The weights rows at chr1:3000 and chr1:4000 describe an A/G SNP.
    The VDS has an A/T SNP at those positions -- the A/G variant does not exist.

    Keying the join on (locus, alleles) means a weights row can only match a
    variant with the same alleles, so both rows are dropped for everyone. This
    is the check the removed path had to opt into via `strict_allele_match`, and
    which it skipped entirely when that was False -- scoring a weight against
    whatever variant happened to occupy the coordinate.
    """
    prs, n_matched = score(vds_lgt, raw_weights, positions=[3000, 4000])

    assert n_matched == 0, "neither A/G variant exists in the VDS"
    assert set(prs.values()) == {0.0}


# --------------------------------------------------------------------------
# Allele orientation. Both of these are Task 2's targets.
# --------------------------------------------------------------------------


def test_silently_drops_rows_of_the_opposite_orientation(vds_lgt, raw_weights):
    """BUG (design limit). `ref_is_effect_allele` is global, not per-row.

    The join key is built from the flag, so a weights file with mixed
    orientation -- some rows naming the REF as the effect allele, some naming
    the ALT -- has whichever half disagrees with the flag silently dropped. It
    is not even counted in `n_matched`. The PGS Catalog does not harmonize the
    effect allele onto the ALT, so mixed files are the normal case here, not a
    pathological one.

    Three of the five rows here vanish with the flag off, four with it on, and
    the score still comes out as a clean number. See TODO.md, Finding 4.
    """
    _, n_alt_oriented = score(vds_lgt, raw_weights, ref_is_effect_allele=False)
    _, n_ref_oriented = score(vds_lgt, raw_weights, ref_is_effect_allele=True)

    assert n_alt_oriented == 2  # 1000, 5000
    assert n_ref_oriented == 1  # 2000 only
    # No setting of the flag scores all five rows.
    assert n_alt_oriented + n_ref_oriented < len(raw_weights.collect())


def test_ref_is_effect_offsets_every_sample_equally(vds_lgt, raw_weights):
    """BUG (benign), and the invariant that makes it benign.

    With `ref_is_effect_allele=True` the code negates the weight and multiplies
    by the ALT count, giving `-w * n_alt` where the truth is `w * (2 - n_alt)`.
    The two differ by a constant `2w` per matched variant. This test pins the
    thing that actually matters: the offset is the *same for every sample*,
    including the hom-ref sample that is never visited. Absolute scores are
    shifted; rankings are not.

    Confirmed on the real All of Us VDS, where every one of 200 samples was
    short by exactly 2w*K. If a change ever makes this offset vary per sample,
    cross-sample comparability is destroyed -- and this test goes red.

    Task 2 removes the offset by adding `2 * sum(w)` over the matched REF-effect
    rows back as a row-level constant. When it does, this test should be
    rewritten to assert the scores equal `truth` outright.
    """
    # chr1:2000 alone. Effect allele is A (the REF). True copies of A:
    # S1 A/A = 2, S2 A/G = 1, S3 G/G = 0.
    truth = {"S1": 2.0, "S2": 1.0, "S3": 0.0}

    prs, n_matched = score(
        vds_lgt, raw_weights, positions=[2000], ref_is_effect_allele=True
    )
    assert n_matched == 1
    assert prs["S1"] == 0.0 and prs["S2"] == -1.0 and prs["S3"] == -2.0

    offsets = {s: truth[s] - prs[s] for s in truth}
    assert set(offsets.values()) == {2.0}, (
        "the 2w offset must be identical for every sample, or scores stop "
        f"being comparable across individuals; got {offsets}"
    )


# --------------------------------------------------------------------------
# Batch must agree with single.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("ref_is_effect_allele", [True, False])
def test_batch_agrees_with_single(vds_lgt, raw_weights, ref_is_effect_allele):
    """The batch calculator gates matches with `hl.if_else` rather than
    `filter_rows`. Different mechanism, same number required."""
    config = PRSConfig(
        ref_is_effect_allele=ref_is_effect_allele, include_n_matched=True
    )
    single, _ = score(
        vds_lgt, raw_weights, ref_is_effect_allele=ref_is_effect_allele
    )

    prepared, _ = _prepare_batch_weights_data({"s1": raw_weights}, config)
    df = _calculate_prs_chunk_batch(
        vds_lgt, {"s1": raw_weights}, prepared, config
    ).to_pandas()
    batch = dict(zip(df["s"], df["s1"], strict=True))

    assert batch == single
