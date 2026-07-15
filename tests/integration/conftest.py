"""Fixtures for the real-hail integration tier.

Nothing here is mocked. A local Spark backend is started, a small GRCh38
VariantDataset is built by hand, and the library's real scoring functions run
against it. `tests/prs/` checks that the code calls the right hail methods;
this tier is the only one that can tell a correct score from a wrong one.

The VDS is built through `hl.vds.VariantDataset.from_merged_representation`
rather than by assembling `reference_data`/`variant_data` directly, because it
reproduces the one property the scoring logic hinges on: a homozygous-reference
sample has **no entry** in `variant_data`. It is not an entry whose genotype is
missing -- it is filtered out of the entry stream entirely, and hail's
aggregators never visit it. See `test_allele_matching.py` for why that matters.
"""

import hail as hl
import pytest

pytestmark = pytest.mark.integration

SAMPLES = ["S1", "S2", "S3", "S4"]

# The mock VDS, one locus per scenario. Genotypes are given as global allele
# indices; the LGT/LA fixture converts them to local indices. A sample absent
# from a row's dict is homozygous reference there and gets a reference block.
#
#   locus         alleles        who carries what
VARIANTS = [
    # Clean biallelic site. Weights name G (the ALT) as the effect allele.
    ("chr1:1000", ["A", "G"], {"S2": [0, 1], "S3": [1, 1]}),
    # Clean biallelic site, but the weights name A (the REF) as the effect
    # allele. S4 is a genuine no-call here: its entry exists, its genotype is
    # missing. That is a different state from S1's absence, and the two are
    # scored differently.
    ("chr1:2000", ["A", "G"], {"S2": [0, 1], "S3": [1, 1], "S4": None}),
    # The VDS has an A/T SNP here. The weights describe an A/G SNP -- a variant
    # that does not exist at this position. Effect allele is A, the REF base.
    ("chr1:3000", ["A", "T"], {"S2": [0, 1], "S3": [1, 1]}),
    # Same mismatch, mirrored: effect allele is G, which is absent from the VDS.
    ("chr1:4000", ["A", "T"], {"S2": [0, 1], "S3": [1, 1]}),
    # Multi-allelic. S2 is C/G, S3 is C/T. Weights name T as the effect allele.
    #
    # This site is also where a REF-effect weight goes wrong if the dosage is
    # taken from the *downcoded* genotype: at the split [C, T] row, S2's G is
    # rewritten to the reference, so its GT is 0/0 and it looks homozygous
    # reference -- but it carries only ONE C, not two. See
    # test_ref_effect_at_a_multiallelic_site.
    ("chr1:5000", ["C", "G", "T"], {"S2": [0, 1], "S3": [0, 2]}),
    # Multi-allelic AND not minimally represented, but the locus does NOT move.
    # A deletion (AGGGC -> A) shares this record with a SNP (AGGGC -> GGGGC).
    # The SNP differs at the *first* base, so `hl.min_rep` trims the shared
    # suffix GGGC and reduces it to A/G at this same locus -- which is how a
    # GWAS names it. Trimming a shared *prefix* would instead move the locus;
    # that case does not occur in All of Us (measured: 0 of 6,001,424 ALTs; see
    # notebooks/measure_minrep_locus_shift.ipynb).
    #
    # S2 carries the SNP, S3 is homozygous for it. No weights row in WEIGHTS
    # names this locus -- the test supplies its own.
    ("chr1:6000", ["AGGGC", "A", "GGGGC"], {"S2": [0, 2], "S3": [2, 2]}),
    # Multi-allelic, and the worst case for a REF-effect weight. S2 is C/C --
    # HOMOZYGOUS for an ALT that the weights row does not name. At the split
    # [A, G] row its genotype downcodes to 0/0, which is indistinguishable from
    # homozygous-reference, so a dosage taken from the downcoded GT would credit
    # S2 with TWO copies of A when it carries NONE.
    ("chr1:7000", ["A", "C", "G"], {"S2": [1, 1], "S3": [0, 2]}),
    # REF sorts AFTER ALT ('G' > 'A'). Every variant above happens to have
    # REF < ALT, which is what let a join keyed on the *sorted* allele pair --
    # against a VDS row key that is NOT sorted -- pass the whole suite while
    # silently dropping ~half of a real weights file. Nothing else here is
    # sensitive to allele order, so this locus is the only thing standing
    # between that bug and a green run. Genotypes mirror chr1:1000 exactly, so
    # the two must score identically under either orientation.
    ("chr1:8000", ["G", "A"], {"S2": [0, 1], "S3": [1, 1]}),
    # BIALLELIC and non-minimally represented: AAAG/GAAG is an A->G SNP written
    # with three shared suffix bases, which is how the real AoU VDS stores some
    # SNPs. `hl.vds.split_multi` min_reps only the rows it SPLITS; an already-
    # biallelic row is passed through with its original alleles. So without an
    # explicit normalization this row stays AAAG/GAAG and never joins a weights
    # row that names it minimally (A/G) -- it scores 0, silently. chr1:6000
    # covers the multi-allelic version, which splitting happens to normalize;
    # this covers the biallelic version, which it does not. Found on the real
    # VDS by validate_scoring_on_aou.ipynb (chr1:1409159).
    ("chr1:8500", ["AAAG", "GAAG"], {"S2": [0, 1], "S3": [1, 1]}),
]

# Every weight is 1.0, so `prs` is literally the summed count of effect-allele
# copies. Score assertions stay readable as copy-number arithmetic.
WEIGHTS = [
    {"chr": "chr1", "pos": 1000, "effect_allele": "G",
     "noneffect_allele": "A", "weight": 1.0},
    {"chr": "chr1", "pos": 2000, "effect_allele": "A",
     "noneffect_allele": "G", "weight": 1.0},
    {"chr": "chr1", "pos": 3000, "effect_allele": "A",
     "noneffect_allele": "G", "weight": 1.0},
    {"chr": "chr1", "pos": 4000, "effect_allele": "G",
     "noneffect_allele": "A", "weight": 1.0},
    {"chr": "chr1", "pos": 5000, "effect_allele": "T",
     "noneffect_allele": "C", "weight": 1.0},
]  # fmt: skip


@pytest.fixture(scope="session", autouse=True)
def hail_context(tmp_path_factory):
    """Starts one local Spark backend for the whole session."""
    hl.init(
        master="local[1]",
        quiet=True,
        skip_logging_configuration=True,
        tmp_dir=str(tmp_path_factory.mktemp("hail")),
    )
    hl.default_reference("GRCh38")
    yield
    hl.stop()


def _build_vds() -> hl.vds.VariantDataset:
    """Builds the mock VDS from VARIANTS.

    Entries use the local (`LGT`/`LA`) encoding, which is what All of Us ships:
    `LA` lists the global allele indices this sample uses, and `LGT` indexes
    into `LA` rather than into `alleles`. `hl.vds.split_multi` converts it to a
    plain `GT` before the scoring code sees it.
    """
    schema = hl.tstruct(
        locus_str=hl.tstr,
        alleles=hl.tarray(hl.tstr),
        s=hl.tstr,
        gt=hl.tarray(hl.tint32),
        la=hl.tarray(hl.tint32),
        GQ=hl.tint32,
        END=hl.tint32,
    )

    entries = []
    for locus_str, alleles, carriers in VARIANTS:
        pos = int(locus_str.split(":")[1])
        for sample in SAMPLES:
            if sample not in carriers:
                # Homozygous reference: a reference block, never a variant
                # entry. This is what makes hom-ref samples invisible to the
                # variant_data entry stream.
                entries.append(
                    {
                        "locus_str": locus_str,
                        "alleles": [alleles[0]],
                        "s": sample,
                        "gt": [0, 0],
                        "la": [0],
                        "GQ": 99,
                        "END": pos,
                    }
                )
                continue

            # LA lists the global allele indices this sample uses; LGT indexes
            # into LA. A no-call keeps LA (and GQ) defined, so the entry still
            # exists -- which is what makes it different from S1's absence.
            global_gt = carriers[sample]
            local_alleles = (
                [0] if global_gt is None else sorted({0, *global_gt})
            )
            local_gt = (
                None
                if global_gt is None
                else [local_alleles.index(i) for i in global_gt]
            )
            entries.append(
                {
                    "locus_str": locus_str,
                    "alleles": alleles,
                    "s": sample,
                    "gt": local_gt,
                    "la": local_alleles,
                    "GQ": 40,
                    "END": None,
                }
            )

    ht = hl.Table.parallelize(entries, schema)
    ht = ht.annotate(
        locus=hl.parse_locus(ht.locus_str, reference_genome="GRCh38")
    )
    ht = ht.annotate(
        LGT=hl.or_missing(hl.is_defined(ht.gt), hl.call(ht.gt[0], ht.gt[1])),
        LA=ht.la,
    )
    ht = ht.drop("gt", "la", "locus_str")

    mt = ht.to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
    vds = hl.vds.VariantDataset.from_merged_representation(mt, is_split=False)
    # Fail loudly on a malformed VDS rather than producing quiet nonsense.
    vds.validate()
    return vds


@pytest.fixture(scope="session")
def vds_lgt():
    """VDS with local (`LGT`/`LA`) genotype encoding -- what AoU ships."""
    return _build_vds()


@pytest.fixture(scope="session")
def raw_weights():
    """The mock GWAS summary, unvalidated -- as a user would hand it over."""
    return hl.Table.parallelize(
        WEIGHTS,
        hl.tstruct(
            chr=hl.tstr,
            pos=hl.tint32,
            effect_allele=hl.tstr,
            noneffect_allele=hl.tstr,
            weight=hl.tfloat64,
        ),
    )


@pytest.fixture(scope="session")
def vds_locus_shifting():
    """A VDS containing a variant whose locus MOVES under `hl.min_rep`.

    `chr1:1001 [GG, G, GT]` -- a deletion of the G at 1002, which VCF anchors
    one base earlier, sharing a record with a G>T SNP at 1002. `min_rep` reduces
    `GG -> GT` to `G -> T` at chr1:**1002**, one base downstream of the row.

    `hl.vds.split_multi` refuses to relocate rows, so it raises on this. That is
    deliberate and this VDS exists to pin it -- see
    `test_a_locus_shifting_variant_raises_rather_than_vanishing`. It is kept out
    of the main fixture because `split_multi` runs over the whole VDS, so a
    single such row would make every other test raise.

    No variant of this shape exists in All of Us: 0 of 6,001,424 ALT alleles
    in a 10Mb window shifted (`notebooks/measure_minrep_locus_shift.ipynb`).
    """
    entries = [
        {
            "locus_str": "chr1:1001",
            "alleles": ["GG", "G", "GT"],
            "s": s,
            "gt": gt,
            "la": [0, 2],
            "GQ": 40,
            "END": None,
        }
        for s, gt in [("S1", [0, 0]), ("S2", [0, 1])]
    ]
    ht = hl.Table.parallelize(
        entries,
        hl.tstruct(
            locus_str=hl.tstr,
            alleles=hl.tarray(hl.tstr),
            s=hl.tstr,
            gt=hl.tarray(hl.tint32),
            la=hl.tarray(hl.tint32),
            GQ=hl.tint32,
            END=hl.tint32,
        ),
    )
    ht = ht.annotate(
        locus=hl.parse_locus(ht.locus_str, reference_genome="GRCh38")
    )
    ht = ht.annotate(LGT=hl.call(ht.gt[0], ht.gt[1]), LA=ht.la)
    ht = ht.drop("gt", "la", "locus_str")
    mt = ht.to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
    return hl.vds.VariantDataset.from_merged_representation(mt, is_split=False)
