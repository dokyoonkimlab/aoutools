"""Indels whose direction the data cannot settle, against a real hail VDS.

At a site whose reference reads `CA...`, the deletion `CA -> C` and the
insertion `CA -> CAA` reduce to `[CA, C]` and `[C, CA]`: the same two alleles in
opposite order, two different variants. A SNP can never be read both ways,
because only one of its alleles can be the reference base; an indel can,
because both of its alleles are prefixes of the reference sequence.

An unordered allele match took a `C`/`CA` weights row onto both split rows. The
deletion row was scored as if its effect allele were the reference, and the
two contributions summed to `2w - w * (copies of the site's other ALTs)`: the
intended insertion canceled out entirely, and a score stayed a clean float.
Found on the real VDS at chr10:100345125 (`[CA, C, CAA, CAAA, TA]`), where
1,152 of a 105,284-variant score's weights were matched twice.

Every weight here is 1.0 (10.0 where two weights must be told apart), so each
expected score is copy-number arithmetic over the genotypes in
`INDEL_PAIR_GENOTYPES` below.
"""

import logging
from collections import Counter

import hail as hl
import pandas as pd
import pytest

from aoutools.prs._calculator import _calculate_prs_chunk, _process_chunks
from aoutools.prs._calculator_batch import (
    _calculate_prs_chunk_batch,
    _prepare_batch_weights_data,
)
from aoutools.prs._calculator_utils import (
    _prepare_weights_for_chunking,
    _validate_and_prepare_weights_table,
)
from aoutools.prs._config import PRSConfig

pytestmark = pytest.mark.integration

# chr1:9000 is the shape found on the real VDS (chr10:100345125 is
# `[CA, C, CAA, CAAA, TA]`). The six samples cover every genotype over its REF
# (global allele 0), deletion (1), and insertion (2). chr1:9100 carries only
# the deletion -- the case a double-match check alone cannot catch.
INDEL_PAIR_GENOTYPES = {
    "S1": None,  # CA/CA -- homozygous reference, no entry
    "S2": [0, 2],  # CA/CAA -- one insertion
    "S3": [2, 2],  # CAA/CAA -- two insertions
    "S4": [0, 1],  # CA/C -- one deletion
    "S5": [1, 2],  # C/CAA -- one of each
    "S6": [1, 1],  # C/C -- two deletions
}
SAMPLES = list(INDEL_PAIR_GENOTYPES)
INDEL_VARIANTS = [
    (
        "chr1:9000",
        ["CA", "C", "CAA"],
        {s: g for s, g in INDEL_PAIR_GENOTYPES.items() if g is not None},
    ),
    ("chr1:9100", ["CA", "C"], {"S4": [0, 1], "S6": [1, 1]}),
]
# Copies of the insertion (global allele 2) and the deletion (allele 1) at
# chr1:9000, per sample.
N_INS = {s: (g or []).count(2) for s, g in INDEL_PAIR_GENOTYPES.items()}
N_DEL = {s: (g or []).count(1) for s, g in INDEL_PAIR_GENOTYPES.items()}
# chr1:9100 carries only the deletion: S4 is CA/C, S6 is C/C.
N_DEL_ONLY_SITE = {s: 0 for s in SAMPLES} | {"S4": 1, "S6": 2}

INSERTION = ("CA", "C")  # effect CA, non-effect C: C -> CA, an insertion
DELETION = ("C", "CA")  # effect C, non-effect CA: CA -> C, a deletion

WEIGHTS_SCHEMA = hl.tstruct(
    chr=hl.tstr,
    pos=hl.tint32,
    effect_allele=hl.tstr,
    noneffect_allele=hl.tstr,
    weight=hl.tfloat64,
)


@pytest.fixture(scope="module")
def vds_indel_pair(build_vds):
    """A VDS with an insertion and a deletion of the same letter at one site,
    and a deletion-only site beside it. See INDEL_VARIANTS."""
    return build_vds(INDEL_VARIANTS, SAMPLES)


def weights(*rows):
    """A weights table from `(pos, (effect, noneffect), weight)` rows."""
    return hl.Table.parallelize(
        [
            {
                "chr": "chr1",
                "pos": pos,
                "effect_allele": effect,
                "noneffect_allele": noneffect,
                "weight": weight,
            }
            for pos, (effect, noneffect), weight in rows
        ],
        WEIGHTS_SCHEMA,
    )


def score(vds, raw_weights, **config_kwargs):
    """Scores through both the single and the batch path; they must agree.

    Returns `({sample: prs}, n_matched, tally)`, the tally being the per-chunk
    indel counts the run's log is built from.
    """
    config = PRSConfig(include_n_matched=True, **config_kwargs)

    tally = Counter()
    prepared = _validate_and_prepare_weights_table(raw_weights, config)
    single = _calculate_prs_chunk(prepared, vds, config, tally=tally)

    batch_tally = Counter()
    grouped, _ = _prepare_batch_weights_data({"s1": raw_weights}, config)
    batch = _calculate_prs_chunk_batch(
        vds, {"s1": raw_weights}, grouped, config, tally=batch_tally
    )

    prs = dict(zip(single["person_id"], single["prs"], strict=True))
    assert dict(zip(batch["person_id"], batch["s1"], strict=True)) == prs
    n_matched = int(single["n_matched"].iloc[0])
    assert int(batch["n_matched_s1"].iloc[0]) == n_matched
    assert {k[1]: v for k, v in batch_tally.items()} == {
        k[1]: v for k, v in tally.items()
    }
    return prs, n_matched, {k[1]: v for k, v in tally.items()}


def test_an_insertion_is_not_also_scored_as_the_deletion(vds_indel_pair):
    """The report's case: an insertion weight at a site that also carries the
    deletion of the same letter scores the insertion, once.

    Under the unordered match this returned `[2, 2, 2, 1, 1, 0]`: every
    insertion carrier scored the same as a homozygous-reference sample, and
    each deletion copy scored -1.
    """
    prs, n_matched, tally = score(vds_indel_pair, weights((9000, INSERTION, 1)))

    assert prs == N_INS
    assert n_matched == 1
    assert tally["n_ambiguous_indel"] == 1
    assert tally["n_reversed_indel"] == 0


def test_a_deletion_is_not_also_scored_as_the_insertion(vds_indel_pair):
    """The mirror case: a deletion weight scores the deletion, once."""
    prs, n_matched, _ = score(vds_indel_pair, weights((9000, DELETION, 1)))

    assert prs == N_DEL
    assert n_matched == 1


def test_insertion_and_deletion_weights_each_score_their_own_variant(
    vds_indel_pair,
):
    """A file that weights both variants at the site keeps both weights.

    Under the unordered match both weights fit both rows, `.find` took the
    first for each row, and the deletion's weight of 10 vanished.
    """
    prs, n_matched, tally = score(
        vds_indel_pair,
        weights((9000, INSERTION, 1.0), (9000, DELETION, 10.0)),
    )

    assert prs == {s: N_INS[s] + 10 * N_DEL[s] for s in SAMPLES}
    assert n_matched == 2
    assert tally["n_ambiguous_indel"] == 2


def test_skip_ambiguous_indels_leaves_them_out(vds_indel_pair):
    """With `skip_ambiguous_indels`, an indel at a site carrying both
    directions contributes nothing and is not counted as matched -- but is
    still counted for the log."""
    prs, n_matched, tally = score(
        vds_indel_pair,
        weights((9000, INSERTION, 1.0), (9000, DELETION, 10.0)),
        skip_ambiguous_indels=True,
    )

    assert prs == {s: 0.0 for s in SAMPLES}
    assert n_matched == 0
    assert tally["n_ambiguous_indel"] == 2


def test_skip_ambiguous_indels_does_not_touch_a_one_direction_site(
    vds_indel_pair,
):
    """Only sites carrying both directions are skipped. The deletion-only site
    scores the same with the option on or off."""
    rows = weights((9100, DELETION, 1.0))

    default, n_default, _ = score(vds_indel_pair, rows)
    skipped, n_skipped, _ = score(
        vds_indel_pair, rows, skip_ambiguous_indels=True
    )

    assert default == skipped == N_DEL_ONLY_SITE
    assert n_default == n_skipped == 1


def test_effect_allele_is_alt_scores_only_the_stated_direction(
    vds_indel_pair,
):
    """`effect_allele_is_alt` states the non-effect allele is the reference,
    so an indel matches only in that order -- and `skip_ambiguous_indels` has
    nothing left to decide."""
    rows = weights((9000, INSERTION, 1.0))

    plain, n_plain, _ = score(vds_indel_pair, rows, effect_allele_is_alt=True)
    with_skip, n_with_skip, _ = score(
        vds_indel_pair,
        rows,
        effect_allele_is_alt=True,
        skip_ambiguous_indels=True,
    )

    assert plain == with_skip == N_INS
    assert n_plain == n_with_skip == 1


def test_effect_allele_is_alt_fast_path_orders_indels(vds_indel_pair):
    """The single-pass path (`effect_allele_is_alt` without `n_matched`) runs
    no row reduction, so no counts and no tripwire -- but the match is the same
    one, and must still score each indel once, in the stated direction."""
    config = PRSConfig(effect_allele_is_alt=True)
    raw = weights((9000, INSERTION, 1.0), (9100, INSERTION, 1.0))

    prepared = _validate_and_prepare_weights_table(raw, config)
    single = _calculate_prs_chunk(prepared, vds_indel_pair, config)
    grouped, _ = _prepare_batch_weights_data({"s1": raw}, config)
    batch = _calculate_prs_chunk_batch(
        vds_indel_pair, {"s1": raw}, grouped, config
    )

    # chr1:9100 has no insertion, so only chr1:9000's insertion scores.
    assert dict(zip(single["person_id"], single["prs"], strict=True)) == N_INS
    assert dict(zip(batch["person_id"], batch["s1"], strict=True)) == N_INS


def test_a_one_direction_site_keeps_the_reference_effect_reading(
    vds_indel_pair,
):
    """Pins a reading the data cannot settle, so a change to it is reviewable.

    chr1:9100 carries only the deletion `[CA, C]`. A weights row with effect
    `CA` and non-effect `C` is an insertion if the file lists the reference as
    the non-effect allele, or this deletion with its reference allele as the
    effect allele if it does not. By default it is kept, as the deletion,
    scored `2 - n_non_ref`: right for a file that lists such indels that way
    (a PGS Catalog file may), and reported as `n_reversed_indel` so a user
    whose file does not can set `effect_allele_is_alt`. A double-match check
    alone would not catch this case, since only one row matches.
    """
    prs, n_matched, tally = score(vds_indel_pair, weights((9100, INSERTION, 1)))

    assert prs == {s: 2 - N_DEL_ONLY_SITE[s] for s in SAMPLES}
    assert n_matched == 1
    assert tally["n_reversed_indel"] == 1
    assert tally["n_ambiguous_indel"] == 0


def test_effect_allele_is_alt_drops_a_reversed_indel(vds_indel_pair):
    """The same row under `effect_allele_is_alt` names an insertion that is not
    in the data, so it matches nothing."""
    prs, n_matched, _ = score(
        vds_indel_pair,
        weights((9100, INSERTION, 1)),
        effect_allele_is_alt=True,
    )

    assert prs == {s: 0.0 for s in SAMPLES}
    assert n_matched == 0


def test_the_counts_are_logged_once_per_run(vds_indel_pair, caplog):
    """The log reports the run's totals once, not once per chunk.

    Two chunks of one weight each: the both-directions insertion at chr1:9000
    (read with the non-effect allele as the reference) and the reversed reading
    at chr1:9100. Both readings are assumptions about the file, so both are
    WARNINGs.
    """
    config = PRSConfig(chunk_size=1)
    full, n_chunks = _prepare_weights_for_chunking(
        weights((9000, INSERTION, 1.0), (9100, INSERTION, 1.0)),
        config,
        validate_table=True,
    )
    assert n_chunks == 2

    with caplog.at_level(logging.INFO, logger="aoutools"):
        partial = _process_chunks(full, n_chunks, vds_indel_pair, config)

    ambiguous = [r for r in caplog.records if "carry both" in r.getMessage()]
    reversed_ = [
        r for r in caplog.records if "effect allele as the" in r.getMessage()
    ]
    assert [(r.levelno, "1 indels" in r.getMessage()) for r in ambiguous] == [
        (logging.WARNING, True)
    ]
    assert [(r.levelno, "1 indels" in r.getMessage()) for r in reversed_] == [
        (logging.WARNING, True)
    ]

    total = pd.concat(partial).groupby("person_id")["prs"].sum().to_dict()
    assert total == {s: N_INS[s] + 2 - N_DEL_ONLY_SITE[s] for s in SAMPLES}


def test_a_weight_matched_to_two_rows_raises(build_vds):
    """The tripwire behind the fix: a weight must never land on two rows.

    Nothing in the current match allows it, so this builds the one shape that
    still could -- the same SNP stored as two rows at one position, once
    minimally (`A/G`) and once padded (`AAAG/GAAG`) -- and checks the run
    raises instead of counting the weight twice.
    """
    vds = build_vds(
        [
            ("chr1:9500", ["A", "G"], {"S1": [0, 1], "S2": [1, 1]}),
            ("chr1:9500", ["AAAG", "GAAG"], {"S1": [0, 1], "S2": [0, 1]}),
        ],
        ["S1", "S2"],
    )
    config = PRSConfig()
    prepared = _validate_and_prepare_weights_table(
        weights((9500, ("G", "A"), 1.0)), config
    )

    with pytest.raises(RuntimeError, match="counted twice"):
        _calculate_prs_chunk(prepared, vds, config)
