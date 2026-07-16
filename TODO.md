# TODO

## Status

`tests/integration/` runs the real scoring code against a real hail VDS
(`pixi run -e integration test-integration`, and a CI job). It exists because the
analysis below was originally done by reading code, and reading code got two
things wrong.

The core findings were then **confirmed on the real All of Us VDS** on the
Workbench. Those numbers are quoted inline below — they are not from the mock.
The live check on current code is `notebooks/validate_scoring_on_aou.ipynb`
(check 1 re-confirms the hom-ref fact); the earlier
`verify_hom_ref_dosage.ipynb`, which measured the now-removed non-split path,
has been retired.

**All three tasks are closed.** The non-split scoring path is gone (Findings 1-3
went with it), allele orientation is resolved per row (Finding 4), and Finding 5
turned out **not to occur in All of Us at all** -- measured, not assumed -- so
Task 3 was cancelled rather than built.

**Finding 6 (fixed).** Found while writing the *oracle* for
`notebooks/validate_scoring_on_aou.ipynb`: at a multi-allelic site, a REF-effect
weight scored a carrier of a *different* ALT as homozygous reference. See below.

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

## ~~Finding 5~~ — the interval prefilter vs. minrep

**Severity was medium-high. CLOSED: the problem does not occur in All of Us.**
Measured with `notebooks/measure_minrep_locus_shift.ipynb`, run on the Workbench.

`hl.min_rep` trims shared bases, and it matters *which end*:

```
minrep of chr1:1001 [AGGGC, A, GGGGC]   (shared SUFFIX)
  AGGGC/GGGGC -> chr1:1001  ['A', 'G']    locus UNCHANGED   <- safe, and relied on

minrep of chr1:1001 [GG, G, GT]         (shared PREFIX)
  GG/GT       -> chr1:1002  ['G', 'T']    locus MOVES +1    <- the feared case
```

Only **prefix** trimming moves a locus, and it happens when a joint caller packs a
SNP into the same record as a deletion whose VCF anchor base sits upstream of it.
The shift equals the distance from that anchor base to the SNP.

The concern was real in principle: `_create_1bp_intervals` builds intervals at the
**weights** locus and `filter_intervals` runs on the **unsplit** VDS, so a row
physically upstream of its minrep'd position would never be read. And
`hl.vds.split_multi` cannot rescue it — hail will not relocate a row, it can only
**raise** (`filter_changed_loci=False`, the default) or **silently drop the
allele** (`True`). Verified both.

**But it does not happen.** In a 10Mb window of chr1 (`chr1:50,000,000-60,000,000`):

| | |
|---|---|
| VDS rows | 4,568,862 |
| multi-allelic rows | 978,586 (**21.4%**) |
| (variant, ALT) pairs | 6,001,424 |
| **ALTs whose locus shifts** | **0** |
| weights rows of PGS000746 recovered by a fix | **0 of 163** |
| weights loci where the crash is reachable | **0** |

All of Us never packs a SNP downstream of an indel's anchor base into one record:
every ALT is already minimally represented against its own row's locus, so
`min_rep` has nothing to relocate. The multi-allelic sites that *do* exist (21% of
rows) are tri-allelic SNPs, or indels sharing an anchor — none of which shift.

### What was done instead of the fix

Nothing was rebuilt. The proposed redesign (replace `hl.vds.split_multi` with an
explode + `min_rep`-keyed join) would recover **zero** variants and would
reintroduce hand-rolled allele handling of exactly the kind that caused Findings
1-3. Cancelled.

Two guards were added instead:

- **`filter_changed_loci` stays at `False` (raise) deliberately.** With a measured
  rate of zero, that exception is not a crash risk — it is a **tripwire**. If a
  future VDS release ever changes variant representation, a shifted variant that
  *reaches* `split_multi` fails loudly instead of quietly dropping and producing a
  plausible, wrong score. It is a **split-step** guard, so it does not cover the
  prefilter case above: a variant whose row sits upstream of its minrep'd,
  GWAS-named locus is dropped by `filter_intervals` *before* the split — a silent
  `n_matched` shortfall, not a raise. So `tests/integration/` and
  `notebooks/validate_synthetic_control_on_aou.ipynb` pin the raise through the
  `_calculate_prs_chunk` seam, which splits the unfiltered VDS. Setting
  `filter_changed_loci` to `True` would convert the tripwire into precisely the
  silent data loss this whole document is about.
  (`test_a_locus_shifting_variant_raises_rather_than_vanishing`)
- **The normalization we *do* rely on is now tested.** Suffix trimming — a
  multi-allelic, non-minimally-represented variant reducing to the `A/G` a GWAS
  actually names — had no coverage at all, despite being the behavior the library
  depends on. (`test_normalizes_a_non_minimal_representation`)

---

## ~~Finding 6~~ — splitting downcodes, and the REF-effect dosage believed it

**Fixed.** This one was live in `dev` after Task 2 and would have shipped.

`hl.vds.split_multi` emits one bi-allelic row per ALT and **downcodes** the
genotype. At the `[C, T]` row of a `C / G,T` site, a sample carrying the G has
its G rewritten to the reference: its `GT` is `0/0` and `GT.n_alt_alleles()` is
`0` — byte-for-byte identical to a genuine homozygous-reference sample.

For an **ALT**-effect weight that is exactly right (the sample really does carry
zero copies of T). For a **REF**-effect weight it is not, because the offset
arithmetic

```
w · (2 − n_alt)  ==  2w − w · n_alt
```

is only valid when `n_alt` counts **every** non-reference allele. After
downcoding it counts one. Measured on the mock VDS with real hail, effect allele
`C` (the REF) at `chr1:5000` (`C / G,T`):

| sample | genotype | true copies of C | scored (before) |
|---|---|---|---|
| S1 | C/C (no entry) | 2 | 2.0 ✅ |
| S2 | **C/G** | **1** | **2.0** ❌ |
| S3 | C/T | 1 | 1.0 ✅ |
| S4 | C/C (no entry) | 2 | 2.0 ✅ |

Worst case is a full `2w`: a sample **homozygous** for the unnamed ALT (`C/C` at
an `A / C,G` site) carries **no** copy of the REF and was credited with two.

**Why it mattered.** 21% of AoU rows are multi-allelic, and the PGS Catalog does
not harmonize the effect allele onto the ALT, so REF-effect rows are the normal
case — not a corner. Every existing test passed, because at a bi-allelic site
`n_alt` and "total non-ref" are the same number.

**The fix.** Extra entry fields survive `split_multi`, so
`_split_multi_with_total_dosage` annotates each entry with `n_non_ref` — the
total non-reference count from its **pre-split** local genotype — before
splitting. `_entry_contribution` counts `n_non_ref` for REF-effect rows and the
downcoded `GT.n_alt_alleles()` for ALT-effect rows. The two branches are *both*
load-bearing; mutating either turns a distinct set of integration tests red and
nothing else.

Pinned by `test_ref_effect_at_a_multiallelic_site`,
`test_ref_effect_when_a_sample_is_homozygous_for_another_alt`, and
`test_alt_effect_is_unaffected_by_another_alt`.

---

## ~~Finding 7~~ — the join key dropped variants two different ways

**Fixed.** Both were silent — a clean float, a low `n_matched`, no error. Both
were live on `dev` after Task 2, shipped in no release, and are now pinned in
`tests/integration/`.

Task 2 keyed the weights on `(locus, hl.sorted([effect, noneffect]))` and joined
against the split `variant_data` row key `(locus, [REF, ALT])`. Two problems:

1. **Allele order.** The VDS row key is *not* sorted, so the sorted weights key
   only met a variant when `REF < ALT`. A `G/A` SNP whose REF sorts after its ALT
   never joined — indistinguishable from a variant absent from the callset. On
   PGS000746 against the real AoU VDS this dropped **922 of 1,940** variants and
   biased the survivors toward `REF ∈ {A, C}`.
2. **Non-minimal representation.** `hl.vds.split_multi` only `min_rep`s the rows
   it *splits*; an already-biallelic row passes through with its original alleles.
   So `AAAG/GAAG` — how the real VDS stores an `A/G` SNP (chr1:1409159) — kept its
   non-minimal alleles and never matched the minimally-named weight. `chr1:6000`
   (multi-allelic, Finding 5) was normalized for free by the split; the biallelic
   passthrough was not.

**The fix — one shuffle-free join.** `_split_multi_with_total_dosage` annotates
every row with `canonical_alleles` (the `min_rep` of its alleles; idempotent on
rows `split_multi` already reduced, so it only changes the biallelic
passthroughs). `_group_weights_by_locus` groups the weights by **locus alone**,
and `_match_weight_at_locus` matches on the **unordered allele set** against
`canonical_alleles`. Keying on locus alone keeps this a key-*prefix* join (no
shuffle); an ordered- or minimal-allele key would force a `key_rows_by` that
shuffles every entry. The set compare is orientation-agnostic, closing (1); the
`canonical_alleles` annotation closes (2). The same `filter_changed_loci`-style
tripwire guards the passthrough annotation (`or_error` on a locus shift), since
`split_multi`'s own guard does not reach the rows it did not split.

Pinned by `test_a_variant_whose_ref_sorts_after_its_alt_scores` (both
orientations) and `test_normalizes_a_non_minimal_biallelic_variant`.

---

## How scoring behaves today

One path, no allele-handling config at all. `hl.vds.split_multi` splits
multi-allelic sites, then the weights are matched to each row **shuffle-free**:
`_group_weights_by_locus` groups them into one array per locus, the key-prefix
join `weights_by_locus[mt.locus]` reads that array with no shuffle, and
`_match_weight_at_locus` picks the row whose **unordered allele set** equals the
row's `canonical_alleles` (its `min_rep`). So the match survives both mismatches
between weights and VDS — allele order and non-minimal representation — with no
re-key of the MatrixTable. See Finding 7 for the two silent drops this replaced.

**Two different dosages**, because splitting downcodes (Finding 6):
`GT.n_alt_alleles()` counts copies of *this row's* ALT; `n_non_ref` counts *all*
non-reference alleles the sample carries at the site, read off the pre-split
genotype.

| weights row | behavior |
|---|---|
| effect allele on the **ALT** | correct: `+w·n_alt`, offset 0; hom-ref contributes nothing, as it should |
| effect allele on the **REF** | correct: `−w·n_non_ref` per entry **plus a row-level `2w`** applied to every sample |
| REF-effect at a multi-allelic site | correct: a carrier of another ALT scores `2 − n_non_ref`, not `2` (Finding 6) |
| variant not in the VDS | dropped for everyone, **and contributes no offset** |
| no-call entry | contributes nothing; the offset is cancelled by `_entry_contribution` |
| non-minimal representation (suffix trim) | correct: normalized to the variant the GWAS names |
| minrep would shift the locus (prefix trim) | **raises** — a deliberate tripwire; does not occur in AoU (Finding 5) |

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
   variant identity *and* disambiguates a file that names two variants at one
   locus. **(Later replaced — the sorted key silently dropped every `REF > ALT`
   variant, ~half of a real file; see Finding 7 for the shuffle-free join that
   supersedes this step.)**
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

### ~~Task 3 — the interval prefilter~~ (cancelled; the problem does not exist)

The prelude measurement came back zero, so the redesign was never built. See
Finding 5 above for the numbers and for the two guards that were added instead.

Worth keeping in mind: this is a property of **All of Us's callset**, not of hail.
If `aoutools` is ever pointed at a VDS from a different pipeline, the assumption
could break — which is the other reason the raising default stays armed.

---

## Workbench validation notebooks (done)

These notebooks close the holes the offline tiers cannot reach. **CI does not run
them** — a human runs them on the Workbench before a release.

- `notebooks/validate_scoring_on_aou.ipynb` — the real-data counterpart of
  `tests/integration/`. Finds a **real variant of each fixture's shape** in the
  real VDS and checks the library against copy numbers counted from
  `hl.vds.to_dense_mt` — an oracle that never calls `min_rep`, `split_multi`, or
  `aoutools`. Writing that oracle is what found **Finding 6**.
- `notebooks/validate_public_api_on_aou.ipynb` — `calculate_prs`,
  `calculate_prs_batch`, and `calculate_pgs` all hard-raise without a `gs://`
  output path, so no offline test reaches them, nor the PGS download, nor
  `_stage_local_file_to_gcs`. This is the only check on any of it.
- `notebooks/validate_synthetic_control_on_aou.ipynb` — the **positive-control**
  tier. **Builds** a synthetic VDS in which every scoring path is present by
  construction, computes the expected PRS independently (pure Python, cross-
  checked against a `to_dense_mt` oracle), and drives the **public API** against
  that known answer — so it checks the user-facing functions and the `gs://`
  round-trip, which the real-PGS notebook cannot verify against ground truth.
- `notebooks/measure_minrep_locus_shift.ipynb` — the measurement behind
  Finding 5: locus-shift rate is **0 of 6,001,424 ALTs** in AoU. Re-run it if the
  split-step tripwire ever fires on a future VDS release.

`aoutools.init_hail()` / `get_vds_path()` (`_workbench.py`) hold the Workbench
boilerplate: requester-pays billing, the reference genome set *after* `hl.init`
(passing `default_reference` to it is deprecated), and a warned fallback for
`$WGS_VDS_PATH`, which current images no longer export. **`DEFAULT_VDS_PATH`
names a specific AoU data release and must be bumped when a new one lands.**

---

## Other open items

- **Release note (Tasks 1 + 2 + Finding 6 — one note, not three).** Every
  previously computed score changes, and users need to know *which kind* of change
  hit them:
  - `split_multi=False` **with REF-effect weights**: the **rankings themselves were
    wrong** (hom-ref samples zeroed, per-sample error). Anything downstream needs
    redoing.
  - `ref_is_effect_allele=True`: rankings were right, absolute values were short by
    a constant. Percentile/z-score results still stand.
  - **Everyone else, including default config**: REF-effect weights rows were
    silently *dropped* and now contribute. Scores move; the old ones were computed
    from a subset of the variants.
  - **Finding 6 — anyone with REF-effect weights at a multi-allelic site**: samples
    carrying a *different* ALT were over-credited by up to `2w` **each**. Per-sample
    and genotype-dependent, so this **reordered the cohort** too. It shipped in no
    release, but it was live on `dev` after Task 2.
  - **Finding 7 — the sorted join key**: every `REF > ALT` variant and every
    non-minimal biallelic variant was silently *dropped* (922 of 1,940 for
    PGS000746), so scores were computed from a biased subset. Also dev-only, live
    after Task 2. Not a concern for released numbers, but noted so the dev-only
    regressions are recorded in one place.

  The "rankings were wrong" bullets are the difference between *your numbers moved*
  and *your result was wrong*. Say it plainly, do not soften it. Three config knobs
  are removed, so this is a breaking change — bump accordingly.
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
