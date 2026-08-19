"""The weights file parser, against real hail.

`tests/prs/test_reader.py` mocks hail entirely, so none of the filtering below
has ever actually run. This tier does, against files on disk.

The case that matters: PGS Catalog harmonized files routinely contain a few
variants that failed to map to the target build, and **every one of them has a
null chromosome AND a null position**. They therefore all collapse onto the same
`variant_id` (`null_null_A_G`) and look like duplicates of each other.

Checking for duplicates before dropping them turned a handful of unmappable rows
into a hard `ValueError` that rejected the entire scoring file -- and because
`calculate_pgs` catches `ValueError` and skips the file, the user lost the whole
score with no visible reason. Real numbers: 20 unmapped of 375,822 in PGS000747,
442 of 3,423,987 in PGS000748, 21 of 1,745,179 in PGS000018. Nearly every large
score has some.
"""

import pytest

from aoutools.prs import _reader
from aoutools.prs._reader import read_prs_weights

pytestmark = pytest.mark.integration

COLUMN_MAP = {
    "chr": "hm_chr",
    "pos": "hm_pos",
    "effect_allele": "effect_allele",
    "noneffect_allele": "other_allele",
    "weight": "effect_weight",
}

HEADER = "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight"


@pytest.fixture(autouse=True)
def no_gcs_staging(monkeypatch):
    """`read_prs_weights` stages local files to GCS, which needs a bucket.

    The parser itself does not care where the file lives, so hand it the local
    path unchanged. This is the only thing mocked in this file.
    """
    monkeypatch.setattr(
        _reader,
        "_stage_local_file_to_gcs",
        lambda file_path, sub_dir: file_path,
    )


def write(tmp_path, rows, name="weights.tsv"):
    path = tmp_path / name
    path.write_text("\n".join([HEADER, *rows]) + "\n")
    return str(path)


def read(path):
    return read_prs_weights(
        file_path=path,
        header=True,
        column_map=COLUMN_MAP,
        delimiter="\t",
        comment="#",
        validate_alleles=True,
        missing="",
        force=True,
    )


def test_unmapped_variants_are_dropped_not_mistaken_for_duplicates(tmp_path):
    """The regression. Five rows: three real, two that failed harmonization.

    The two unmapped rows share a `variant_id` of `null_null_...`, purely
    because both of their coordinates are null. They are not duplicates of each
    other in any meaningful sense, and they must not abort the file -- they must
    be dropped, leaving the three scoreable variants.
    """
    path = write(
        tmp_path,
        [
            "1\t1000\tA\tG\t0.1",
            "1\t2000\tC\tT\t0.2",
            "2\t3000\tG\tA\t0.3",
            "\t\tA\tG\t0.4",  # failed to map
            "\t\tC\tT\t0.5",  # failed to map; same null_null id shape
        ],
    )

    ht = read(path)

    assert ht.count() == 3, (
        "the three mappable variants must survive; the two unmapped rows must "
        "be dropped rather than rejecting the whole file as duplicated"
    )
    positions = sorted(ht.pos.collect())
    assert positions == [1000, 2000, 3000]


def test_a_genuine_duplicate_still_raises(tmp_path):
    """The duplicate check must still work for real duplicates.

    The fix drops unmapped rows *before* the check. It must not weaken the
    check itself: a variant named twice at the same locus with the same alleles
    would be double-counted, so it is still an error.
    """
    path = write(
        tmp_path,
        [
            "1\t1000\tA\tG\t0.1",
            "1\t1000\tA\tG\t0.2",  # same locus, same alleles
            "1\t2000\tC\tT\t0.3",
        ],
    )

    with pytest.raises(ValueError, match="[Dd]uplicate"):
        read(path)


def test_all_variants_unmapped_is_an_error(tmp_path):
    """A file with nothing scoreable in it must not return an empty table.

    Silently handing back zero variants would produce a PRS of 0.0 for every
    sample -- a clean, plausible, meaningless column.
    """
    path = write(tmp_path, ["\t\tA\tG\t0.1", "\t\tC\tT\t0.2"])

    with pytest.raises(ValueError, match="coordinates|empty|filtered out"):
        read(path)


def test_missing_and_zero_weights_are_still_dropped(tmp_path):
    """Pre-existing behavior, now covered by a real parse."""
    path = write(
        tmp_path,
        [
            "1\t1000\tA\tG\t0.1",
            "1\t2000\tC\tT\t0",  # zero weight: contributes nothing
            "1\t3000\tG\tA\t",  # missing weight
        ],
    )

    ht = read(path)

    assert ht.count() == 1
    assert ht.pos.collect() == [1000]


def test_chromosome_prefix_is_standardized(tmp_path):
    """The VDS uses `chr1`; weights files often say `1`. Mixed input is fine."""
    path = write(
        tmp_path,
        [
            "1\t1000\tA\tG\t0.1",
            "chr2\t2000\tC\tT\t0.2",
        ],
    )

    ht = read(path)

    assert sorted(ht.chr.collect()) == ["chr1", "chr2"]


def test_lowercase_alleles_are_uppercased_on_import(tmp_path):
    """Allele case is normalized at the door, not treated as invalid data.

    Everything downstream compares allele strings literally -- the ACGT check
    in `_validate_alleles`, `_check_duplicated_ids`, and the VDS join -- and
    hail compares strings case-sensitively. A lowercase file used to survive
    import untouched and then match nothing in the VDS, scoring every sample
    0.0 in silence. See
    `test_allele_matching.py::test_lowercase_alleles_score_the_same_as_uppercase`.

    Note this runs with `validate_alleles=True` (the `read` helper above). The
    normalization deliberately happens *before* that check, so a lowercase row
    is accepted rather than dropped against `^[ACGT]+$` -- 'a' is a valid base,
    just written in a different case.
    """
    path = write(
        tmp_path,
        [
            "1\t1000\ta\tg\t0.1",
            "1\t2000\tC\tt\t0.2",  # mixed case within one row
            "1\t3000\tG\tA\t0.3",  # already uppercase, must pass through
        ],
    )
    table = read(path).collect()

    assert len(table) == 3, "no row may be dropped for being lowercase"
    assert [(r.effect_allele, r.noneffect_allele) for r in table] == [
        ("A", "G"),
        ("C", "T"),
        ("G", "A"),
    ]


def test_alleles_differing_only_in_case_are_caught_as_duplicates(tmp_path):
    """The reason the reader normalizes *before* `_check_duplicated_ids`.

    `chr1:1000 a/g` and `chr1:1000 A/G` are the same variant written twice. If
    case were normalized later -- in `_validate_and_prepare_weights_table`, on
    the calculator side -- the duplicate check would compare the raw strings,
    see two distinct `variant_id`s, and pass. The scorer would then collect both
    rows at the locus and `_match_weight_at_locus`'s `.find` would silently use
    only the first, dropping a weight with no warning.

    Normalizing at import turns that into the loud `ValueError` it should be.
    """
    path = write(
        tmp_path,
        [
            "1\t1000\ta\tg\t0.1",
            "1\t1000\tA\tG\t0.2",
        ],
    )
    with pytest.raises(ValueError, match="Duplicate variants found"):
        read(path)
