# TODO

## Status

`tests/integration/` runs the real scoring code against a real hail VDS
(`pixi run -e integration test-integration`, and a CI job). It exists because the
analysis below was originally done by reading code, and reading code got two
things wrong.

The core findings were then **confirmed on the real All of Us VDS** with
`notebooks/verify_hom_ref_dosage.ipynb`, run on the Workbench. Those numbers are
quoted inline below — they are not from the mock.

**Tasks 1 and 2 are done.** The non-split scoring path is gone (Findings 1-3 went
with it) and allele orientation is now resolved per row (Finding 4). **Task 3 is
open**, and Finding 5 is still live.

---

## The root cause: a hom-ref sample is not a missing entry

Everything else follows from this. In a VDS, `variant_data` holds entries **only**
for samples with a non-reference call at a row. A homozygous-reference sample is
absent — and hail *filters absent entries out of the entry stream*. Aggregators
never visit them.

**Confirmed on All of Us** (5 real PCSK9 SNPs, 200 samples): 94 samples have **zero
entries** in `variant_data` yet all 5 entries in the dense matrix. Across the window
there are 113 entries where a dense matrix would hold 1,000 — ~89% of the genotype
stream is simply absent. This is the normal state of a VDS, not an edge case.
(`test_hom_ref_samples_are_absent_from_the_entry_stream` pins it in the mock.)

The practical consequence, and the thing to keep in mind for Task 2: **you cannot
give a hom-ref sample a dosage by handling a missing genotype.** It has no entry to
handle. Anything a hom-ref sample must contribute has to arrive as a *row-level
constant*, applied to every sample regardless of genotype.

---

## ~~Finding 4~~ — `ref_is_effect_allele` was a file-level flag for a row-level property

**Severity was high. FIXED by Task 2.** Confirmed against the PGS Catalog format spec.

The join key is built from a **global** flag, but allele orientation is a **per-row**
property. Every weights row whose orientation disagrees with the flag produces a key
that does not exist in the VDS and is **silently dropped** — not scored, and not even
counted in `n_matched`.

The assumption that PGS Catalog files are harmonized so that effect = ALT is
**false**. Their format spec says of `effect_allele`:

> "this does not necessarily need to correspond to the minor allele/alternative
> allele"

and of `other_allele`: *"this does not necessarily need to correspond to the
reference allele."* The harmonized (HmPOS) files remap coordinates, rsIDs, and
strand — they explicitly **preserve** the original effect/other allele columns. So
"harmonized to GRCh38" guarantees the coordinates, not the orientation.

This is expected, not pathological: the effect allele is whichever allele the effect
size was estimated for (usually risk-increasing or minor), while REF/ALT is a
property of the reference assembly. The two disagree wherever the reference carries
the minor allele. On the 5-row mock file, the flag off matches 2 rows and on matches
1 — **no setting scores them all**, and the score still comes out a clean float.
(`test_silently_drops_rows_of_the_opposite_orientation`)

Related, same fix: with `ref_is_effect_allele=True` the code computed `-w·n_alt`
where the truth is `w·(2 - n_alt)`, so every score was short by a constant `2w` per
matched variant. The Workbench run gave `10.0` for all 200 samples (`2w·K`, `w=1`,
`K=5`) — a single constant, so rankings survived but absolute scores were
meaningless.

**Both are fixed.** See Task 2 below for what replaced them.

---

## Finding 5 — the 1bp interval prefilter defeats the split path's own minrep

**Severity: medium-high. Live; no test covers it yet.**

`hl.vds.split_multi` applies `hl.min_rep`, which can **move the locus**. Verified on
real hail:

```
minrep of chr1:1000 [AGG, AG, AGT]:
  AGG/AG   -> chr1:1000  ['AG', 'A']    locus unchanged
  AGG/AGT  -> chr1:1002  ['G',  'T']    locus MOVES +2
```

A GWAS reports that SNP at its normalized position, chr1:**1002**. But
`_create_1bp_intervals` (`_calculator_utils.py`) builds the interval at the
**weights** locus, and `hl.vds.filter_intervals` runs on the **unsplit** VDS
(`_calculator.py`) — *before* `split_multi`, where the row still sits at chr1:1000:

```
1bp interval [chr1:1002] -> VDS rows kept: 0    <-- the row lives at 1000
1bp interval [chr1:1000] -> VDS rows kept: 1
```

The variant is discarded before `split_multi` can normalize it.

Note the limit of this: minrep only shifts the locus when the two alleles share a
**leading** base. VDS alleles `['AGGGC', 'A', 'GGGGC']` minrep to `['A', 'G']` at
the *same* locus, and that variant matches fine. Only the shared-prefix case is lost.

Two faces:

- **Silent miss** (the common one): weights locus ≠ VDS locus → row filtered out →
  variant never scores, no error.
- **Hard crash** (latent): `hl.vds.split_multi(vds)` is called with the default
  `filter_changed_loci=False`, which **raises** rather than filters when a REF/ALT
  pair changes locus. If a weights locus lands exactly on the *start* of a padded
  multi-allelic, the row survives the interval filter, `split_multi` runs on it, and
  the run dies. Never observed in practice — because the prefilter usually throws
  those rows away first. **Widening the intervals (the fix) makes this reachable**,
  so `filter_changed_loci` must be settled in the same change.

---

## How scoring behaves today

One path, no allele-handling config at all. `hl.vds.split_multi` splits
multi-allelic sites, the weights are joined on `(locus, **sorted** allele pair)` —
so the join key *is* the allele check, and it matches either orientation — and
dosage is `GT.n_alt_alleles()`.

| weights row | behavior |
|---|---|
| effect allele on the **ALT** | correct: `+w·n_alt`, offset 0; hom-ref contributes nothing, as it should |
| effect allele on the **REF** | correct: `−w·n_alt` per entry **plus a row-level `2w`** applied to every sample |
| variant not in the VDS | dropped for everyone, **and contributes no offset** |
| no-call entry | contributes nothing; the offset is cancelled by `_entry_contribution` |
| minrep shifts the locus | **silently missed** (Finding 5, still open) |

---

## The plan

### ~~Task 1 — remove the non-split path~~ (done)

`split_multi` and `strict_allele_match` are gone from `PRSConfig`, along with
`_prepare_mt_non_split`, `_calculate_dosage`, `_check_allele_match`, and the batch
gates. Passing either knob is now a `TypeError`
(`tests/prs/test_config.py::test_removed_non_split_params_are_rejected`).

This removed three findings outright:

- **Finding 1 (high, ranking-corrupting).** With the effect allele on the REF base,
  the non-split score degenerated to *"how many heterozygous sites does this person
  have"* — a hom-ref sample carries 2 copies and was never visited, a hom-alt carries
  0 and also scored 0. On the Workbench run the error was exactly `2 × n_hom_ref`,
  ranging over `[4, 6, 8, 10]` across samples: **per-sample and genotype-dependent,
  not an offset.** Sample `1000004` truly carried 8 copies and scored **0.0**; sample
  `1000774`, carrying 7, scored **3.0** and ranked above them. Spearman against truth
  was < 1.
- **Finding 2 (low, theoretical).** A no-call entry — genotype present but unknown —
  was handed 2 copies of the effect allele. The Workbench run found 0 no-call entries
  out of 113, so this was never biting in practice. Now impossible:
  `test_does_not_invent_a_genotype_for_a_no_call`.
- **Finding 3 (medium).** `strict_allele_match=False` did no allele check at all, so
  a weights row was scored against whatever variant occupied its coordinate — the
  effect allele matching the REF base is satisfied by *any* variant at that locus.

Removal rather than repair, because the non-split path was not merely "faster but
less robust": it did **no normalization at all**, so a GWAS `G/T` row could never
match VDS alleles `[AGG, AG, AGT]` at any interval width. It was incapable, not just
slower. Its speed advantage was never measured.

**Note:** removing it did **not** gain the multi-allelic matching described in
Finding 5 — the split path drops those variants too, for a different reason. That
arrives only with Task 3.

### ~~Task 2 — fix allele orientation~~ (done)

`ref_is_effect_allele` is gone. Orientation is resolved **per row**, against the
VDS's own REF/ALT:

1. `_key_weights_by_variant` keys the weights on `(locus, hl.sorted([effect,
   noneffect]))`. Sorting canonicalizes the pair, so a row matches its variant
   whichever way round it is written; keeping the alleles in the key preserves
   variant identity *and* disambiguates a file that names two variants at one locus.
2. `_orient_weight_and_offset` returns `(weight_per_alt_copy, hom_ref_offset)` —
   `(+w, 0)` when the effect allele is the ALT, `(−w, 2w)` when it is the REF.
3. `prs = agg.sum(contribution) + agg_rows.sum(hom_ref_offset)`.

The algebra: a REF-effect row's true contribution is `w·(2 − n_alt) = 2w − w·n_alt`.
The `−w·n_alt` term is **already correct for absent samples** (`n_alt = 0` → 0), so
the aggregator skipping them costs nothing. Only the `2w` is missing, and it does not
depend on genotype — it is a row-level scalar. `mt.aggregate_rows(..., _localize=False)`
keeps it a lazy expression, so it folds into the existing job: **no densification, no
second entry pass, no extra VDS read.**

Three invariants hold this up, and all three are asserted:

- **The offset only applies to matched rows.** Crediting `2w` for a variant that is
  not in the callset would invent signal for the whole cohort. The single-score path
  gets this from `filter_rows`; the batch path never filters rows, so it gates on
  `is_valid_{score}`. (`test_does_not_credit_the_offset_for_an_unmatched_variant`,
  `test_batch_agrees_with_single`)
- **The offset is cancelled for a no-call.** A hom-ref sample has no entry; a no-call
  *does*, so it would collect the `2w` too — reintroducing Finding 2 by a new route.
  `_entry_contribution` gives a no-call `−hom_ref_offset`, zeroing it out. Without
  this, whether an unknown genotype was dropped or imputed would depend on which
  allele the GWAS happened to label the effect allele, which is arbitrary.
  (`test_does_not_invent_a_genotype_for_a_no_call`)
- **It survives chunking**, summing once per matched row across chunks.
  (`test_chunking_partitions_the_weights`)

Both were **mutation-tested**: dropping the offset, and dropping the no-call
cancellation, each turn exactly the expected tests red and nothing else.

### Task 3 — the interval prefilter (Finding 5)

**Prelude (do this first): confirm on the real AoU VDS, in a notebook.** Everything
below is sized by the answer, and the mock cannot tell us the rate. Measure:

1. How often does a variant in the AoU `variant_data` have a REF longer than 1bp
   *and* a multi-allelic ALT whose minrep shifts the locus? (Scan a window; count.)
2. For a real PGS Catalog weights file: how many loci fail to match the VDS today,
   and how many of those are explained by a locus-shifted multi-allelic sitting
   upstream within N bp? Sweep N to pick a pad width empirically.
3. Does any weights locus land exactly on the start of such a multi-allelic — i.e.
   is the `filter_changed_loci=False` crash reachable on real data? (The user has
   never hit it, which is consistent with the prefilter hiding it.)

Then implement: left-pad the intervals in `_create_1bp_intervals` by the measured
width so the padded multi-allelic survives to `split_multi`, and settle
`filter_changed_loci` (widening the intervals makes the crash *more* reachable, not
less — a variant that shifts out of its interval must be handled deliberately, not by
an exception). Padding costs VDS read, which is real money on this dataset, so the
width must be justified by (2), not guessed.

---

## Other open items

- **Release note (Tasks 1 + 2 together — one note, not two).** Every previously
  computed score changes, and users need to know *which kind* of change hit them:
  - `split_multi=False` **with REF-effect weights**: the **rankings themselves were
    wrong** (hom-ref samples zeroed, per-sample error). Anything downstream needs
    redoing.
  - `ref_is_effect_allele=True`: rankings were right, absolute values were short by
    a constant. Percentile/z-score results still stand.
  - **Everyone else, including default config**: REF-effect weights rows were
    silently *dropped* and now contribute. Scores move; the old ones were computed
    from a subset of the variants.

  That first bullet is the difference between "your numbers moved" and "your result
  was wrong". Say it plainly, do not soften it. Three config knobs are removed, so
  this is a breaking change — bump accordingly.
- **`main` is behind `dev`.** All current work is on `dev`. Clean fast-forward
  whenever wanted.
- **No strand harmonization exists anywhere.** Weights are assumed to be on the same
  strand and build as the VDS. Palindromic variants (A/T, C/G) are undetected. Out of
  scope for Tasks 2-3, but it is the most likely cause of a non-matching weights row.
- **`_standardize_chromosome_column` and scaffold contigs.** The vectorized rewrite
  (`0c6f467`) prefixes every unprefixed contig, so GRCh38 scaffold/alt contigs
  (`KI270728.1`) become `chrKI270728.1` and make `hl.locus` throw. Fails loudly, and
  PGS Catalog / PRS-CS files don't emit scaffold contigs, so this is theoretical —
  but it is a narrow regression.
- **Bundled weights data is still unused by tests.** `aoutools/data/*.{csv,tsv}` are
  100 real GRCh38 chr1 variants whose `ID` column encodes `chr:pos:REF:ALT`. They
  would give `_reader.py` its first real-hail coverage (100% mocked today) and
  exercise chunking at a more realistic scale. Note they are uniformly
  `A1 = effect = ALT`, which is a PRS-CS output convention — not a property of
  weights files in general (Finding 4).
