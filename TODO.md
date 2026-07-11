# TODO

## Status

`tests/integration/` runs the real scoring code against a real hail VDS
(`pixi run -e integration test-integration`, and a CI job). It exists because the
analysis below was originally done by reading code, and reading code got two
things wrong.

The core findings were then **confirmed on the real All of Us VDS** with
`notebooks/verify_hom_ref_dosage.ipynb`, run on the Workbench. Those numbers are
quoted inline below — they are not from the mock.

**No `aoutools/` code has been changed.** The findings are characterized, not
fixed. The three tasks at the bottom are the agreed plan.

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

`_calculate_dosage` (`_calculator_utils.py:214`) opens with a branch that treats a
missing genotype as homozygous reference — dosage `2` if the effect allele is REF —
documented as *"accounting for sparse storage of homozygous reference calls"*.

**That branch never sees a sparse hom-ref call.** It cannot: those entries are not
in the stream. It fires for exactly one thing — an entry that *exists* with a
missing genotype, i.e. a genuine **no-call**. It is precisely inverted from its
stated purpose.

---

## Finding 1 — non-split + effect allele on the REF loses every hom-ref sample

**Severity: high. The only finding that moves rankings, not just absolute scores.**
**Confirmed on All of Us.** Config: `split_multi=False` (any `strict_allele_match`).

With a REF effect allele the score degenerates to *"how many heterozygous sites does
this person have"*: a hom-ref sample carries 2 copies and is never visited (scores
0), a hom-alt carries 0 and also scores 0.

On the Workbench run, the error was exactly `2 × n_hom_ref` for every sample,
ranging over `[4, 6, 8, 10]` — **per-sample and genotype-dependent, not an offset**.
Sample `1000004` truly carries 8 copies and scored **0.0**; sample `1000774`,
carrying 7, scored **3.0** and ranked above them. Spearman against truth < 1.

Contrast the split path with `ref_is_effect_allele=True`: it computes `-w·n_alt`
instead of `w·(2 - n_alt)`, off by `2w` per matched variant — but off by `2w` for
**every sample including the hom-ref one**. The Workbench run gave `10.0` for all
200 samples (`2w·K`, `w=1`, `K=5`): a single constant. Absolute scores shift,
rankings survive. (`test_split_ref_is_effect_offsets_every_sample_equally`)

So: the split path's offset is benign; the non-split path's is not.

---

## Finding 2 — a no-call is scored as two copies of the reference

**Severity: low — theoretical.** A `variant_data` entry that exists with a missing
genotype is visited, the missing-genotype branch fires, and the sample is handed
**2 copies** of the effect allele. An unknown genotype should be excluded, not
imputed. (`test_non_split_scores_a_no_call_as_hom_ref`)

**The Workbench run found 0 no-call entries out of 113.** AoU's `variant_data` does
not appear to carry them, at least in that window. Real, but not currently biting.

---

## Finding 3 — `strict_allele_match=False` scores variants the GWAS never studied

**Severity: medium.** Requires opting out of two safe defaults (`_config.py:70-72`).

The dosage arithmetic is **not** wrong — `_calculate_dosage` string-matches the
effect allele truthfully. What is lost is **variant identity**. With
`strict_allele_match=False` there is no allele check at all; the only filter is
`hl.is_defined(mt.weights_info)` on a locus-only join. And `effect_allele ==
alleles[0]` is satisfied by *any* variant at that locus, because the reference base
is the same string whichever ALT is there.

Weights say `chr1:3000`, effect `A`, noneffect `G` — an A/G SNP. The VDS has an
**A/T** SNP there; the A/G SNP does not exist. A sample who is A/T picks up one
copy of `A`, purely because of an unrelated `T` genotype. The weight was estimated
for an A-vs-G contrast. Per-sample noise, so it perturbs rankings.
(`test_non_split_loose_scores_a_variant_the_gwas_never_studied`)

The mirror case fails safe: effect allele is an ALT that isn't present → dosage 0 →
no contribution, though the row still counts in `n_matched`.

If a GWAS names an A/G SNP and AoU has no such variant, the likely cause is a **bad
weights row** — strand flip (no strand harmonization exists anywhere; palindromic
A/T and C/G are undetected), build mismatch, bad rsID→allele mapping — not a gap in
a callset that size. Dropping the row, which strict does, is the right action.

---

## Finding 4 — `ref_is_effect_allele` is a file-level flag for a row-level property

**Severity: high. Confirmed against the PGS Catalog format spec.**

The split path builds its `(locus, alleles)` join key from a **global** flag, but
allele orientation is a **per-row** property. Every weights row whose orientation
disagrees with the flag produces a key that does not exist in the VDS and is
**silently dropped** — not scored, and not even counted in `n_matched`.

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
(`test_split_silently_drops_rows_of_the_opposite_orientation`)

---

## Finding 5 — the 1bp interval prefilter defeats the split path's own minrep

**Severity: medium-high. Newly found; no test covers it yet.**

`hl.vds.split_multi` applies `hl.min_rep`, which can **move the locus**. Verified on
real hail:

```
minrep of chr1:1000 [AGG, AG, AGT]:
  AGG/AG   -> chr1:1000  ['AG', 'A']    locus unchanged
  AGG/AGT  -> chr1:1002  ['G',  'T']    locus MOVES +2
```

A GWAS reports that SNP at its normalized position, chr1:**1002**. But
`_create_1bp_intervals` (`_calculator_utils.py:355`) builds the interval at the
**weights** locus, and `hl.vds.filter_intervals` runs on the **unsplit** VDS
(`_calculator.py:253-258`) — *before* `split_multi` (`_calculator.py:62`), where the
row still sits at chr1:1000:

```
1bp interval [chr1:1002] -> VDS rows kept: 0    <-- the row lives at 1000
1bp interval [chr1:1000] -> VDS rows kept: 1
```

The variant is discarded before `split_multi` can normalize it. `_calculator.py:336-344`
explicitly claims this matching works. It does not.

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

## How each config counts the effect allele (current behavior)

| config | effect = **ALT** rows | effect = **REF** rows |
|---|---|---|
| `split=True, ref_is_effect=False` (default) | correct: `+w·n_alt` | **silently dropped** (Finding 4) |
| `split=True, ref_is_effect=True` | **silently dropped** (Finding 4) | `-w·n_alt`; flat `Σ2w` offset, rankings safe |
| `split=False, strict=True` | correct | **hom-ref samples lost → reorders people** (Finding 1) |
| `split=False, strict=False` | correct if the variant is real; 0 if not (fails safe) | same reordering, **plus** phantom variants (Finding 3) |

In every config, a sample absent from `variant_data` contributes nothing, and any
variant whose minrep shifts its locus is dropped before scoring (Finding 5).
`n_matched` counts rows that survived the row filter, so under `strict=False` a row
contributing dosage 0 still counts as "matched" — it over-reports overlap.

---

## The plan: three tasks

### Task 1 — remove the non-split path

Delete `split_multi` and `strict_allele_match` from `PRSConfig`, and with them
`_prepare_mt_non_split`, `_calculate_dosage` (both `GT` and `LGT`/`LA` branches),
`_check_allele_match`, and the batch gates (`_calculator_batch.py:141,194,246`).

**Kills Findings 1, 2, and 3 outright**, and removes a second dosage implementation
that has to be kept in sync with the first.

Why removal rather than repair — the non-split path is not merely "faster but less
robust":

- It does **no normalization at all**. String-matching a GWAS `G/T` row against VDS
  alleles `[AGG, AG, AGT]` can never match, at any interval width. On this class of
  variant it is not slower-but-equivalent; it is incapable.
- It owns the only bug that reorders a cohort.
- The speed advantage is undocumented and unmeasured. `hl.vds.split_multi` is a
  per-row explode, not a shuffle.

Rewrite the integration tests that currently pin Findings 1/2/3 as removal proofs.
Self-contained; unblocked.

**Note:** removing this path does **not** by itself gain the multi-allelic matching
described above — Finding 5 means the split path drops those variants too, for a
different reason. That gain arrives only with Task 3.

### Task 2 — fix allele orientation on the split path

Remove `ref_is_effect_allele` entirely. It is a global flag standing in for a
per-row property (Finding 4).

Post-`split_multi` rows are biallelic, so:

1. Key the weights on **locus** and identify the variant orientation-free:
   `hl.set(alleles) == hl.set([effect_allele, noneffect_allele])`.
2. Per row: if `effect == alleles[1]`, contribute `+w·n_alt`. Otherwise (effect is
   REF) contribute `−w·n_alt` **and add `2w` to a row-level constant.**
3. `prs = agg.sum(w_eff · n_alt) + 2·Σw` over matched REF-effect rows.

The algebra: the true contribution of a REF-effect row is `w·(2 − n_alt) = 2w −
w·n_alt`. The `−w·n_alt` term is **already correct for absent samples** (`n_alt = 0`
→ 0), so the aggregator skipping them costs nothing. Only the `2w` is missing, and
it does not depend on genotype — it is a row-level scalar.

The correction is a **rows-only aggregation** over already-filtered rows: no
densification, no second entry pass, no extra VDS read. This supersedes the earlier
belief that fixing the offset required an expensive densify.

Result: mixed-orientation weights files score every row, the `2w` offset is gone, and
two config knobs disappear.

### Task 3 — Finding 5: the interval prefilter

**Prelude (do this first): confirm on the real AoU VDS, in a notebook.** Everything
below is sized by the answer, and the mock cannot tell us the rate. Measure:

1. How often does a variant in the AoU `variant_data` have a REF longer than 1bp
   *and* a multi-allelic ALT whose minrep shifts the locus? (Scan a window; count.)
2. For a real PGS Catalog weights file: how many loci fail to match the VDS today,
   and how many of those are explained by a locus-shifted multi-allelic sitting
   upstream within N bp? Sweep N to pick a pad width empirically.
3. Does any weights locus land exactly on the start of such a multi-allelic — i.e.
   is the `filter_changed_loci=False` crash reachable on real data?

Then implement: left-pad the intervals in `_create_1bp_intervals` by the measured
width so the padded multi-allelic survives to `split_multi`, and settle
`filter_changed_loci` (widening the intervals makes the crash *more* reachable, not
less — a variant that shifts out of its interval must be handled deliberately, not
by an exception). Padding costs VDS read, which is real money on this dataset, so
the width must be justified by (2), not guessed.

Blocked on nothing, but do it after Tasks 1–2 so it lands on the simplified code.

---

## Other open items

- **`main` is behind `dev`.** All current work is on `dev`. Clean fast-forward
  whenever wanted.
- **No strand harmonization exists anywhere.** Weights are assumed to be on the same
  strand and build as the VDS. Palindromic variants (A/T, C/G) are undetected. Out
  of scope for the three tasks, but it is the most likely cause of a non-matching
  weights row.
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
