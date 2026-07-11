# TODO

## Status

`tests/integration/` now runs the real scoring code against a real hail VDS
(`pixi run -e integration test-integration`, and a CI job). It exists because the
analysis below was originally done by reading code, and reading code got two
things wrong. Every claim here is now pinned by a passing test.

**No `aoutools/` code has been changed.** The findings are characterized, not
fixed. The decisions at the bottom are open.

---

## The root cause: a hom-ref sample is not a missing entry

Everything else follows from this. In a VDS, `variant_data` holds entries **only**
for samples with a non-reference call at a row. A homozygous-reference sample is
absent — and hail *filters absent entries out of the entry stream*. Aggregators
never visit them.

```
agg.count() per sample over variant_data, 5 rows, 4 samples:
  S1 (hom-ref everywhere)  -> 0      <-- not "visited and missing". Not visited.
  S2 (a call at each site) -> 5
  S3 (a call at each site) -> 5
  S4 (one no-call entry)   -> 1
```
(`test_hom_ref_samples_are_absent_from_the_entry_stream`)

`_calculate_dosage` (`_calculator_utils.py:214`) opens with a branch that treats a
missing genotype as homozygous reference — dosage `2` if the effect allele is REF —
documented as *"accounting for sparse storage of homozygous reference calls"*.

**That branch never sees a sparse hom-ref call.** It cannot: those entries are not
in the stream. The branch fires for exactly one thing — an entry that *exists* with
a missing genotype, i.e. a genuine **no-call**. It is precisely inverted from its
stated purpose, and both halves of the inversion are bugs.

---

## Finding 1 — non-split + effect allele on the REF loses every hom-ref sample

**Severity: high. This one moves rankings, not just absolute scores.**

Config: `split_multi=False` (any `strict_allele_match`). Weights row whose effect
allele is the reference base — common in real PGS Catalog files.

At `chr1:2000`, VDS alleles `A/G`, effect allele `A`:

| sample | genotype | true copies of A | scored | |
|---|---|---|---|---|
| S1 | A/A | **2** | **0** | never visited |
| S2 | A/G | 1 | 1 | correct |
| S3 | G/G | 0 | 0 | correct |

The true ordering `S1 > S2 > S3` comes out as `S2 > S1 = S3`. The error applies
only to samples who are hom-ref at the site, so it is **genotype-dependent and
per-sample** — not a uniform offset. Two people are scored differently for reasons
unrelated to their genotype at the variant.

Contrast the split path, which gets this right in the way that matters: with
`ref_is_effect_allele=True` it computes `-w * n_alt` instead of `w * (2 - n_alt)`,
which is off by a constant `2w` — but off by `2w` for **every sample including the
hom-ref one**, because the hom-ref sample contributes 0 to both. Absolute scores
shift; rankings survive. That invariant is now asserted directly
(`test_split_ref_is_effect_offsets_every_sample_equally`).

So: the split path's offset is benign; the non-split path's is not.

---

## Finding 2 — a no-call is scored as two copies of the reference

**Severity: medium. Applies on the default config too.**

S4's entry at `chr1:2000` exists with a missing genotype — an *unknown* genotype.
It is visited, the missing-genotype branch fires, and S4 is handed **2 copies** of
the effect allele. An unknown genotype should be excluded, not imputed to hom-ref.
(`test_non_split_scores_a_no_call_as_hom_ref`)

Worth confirming against a real AoU VDS whether no-call entries actually occur in
`variant_data`. If they don't, this is theoretical.

---

## Finding 3 — `strict_allele_match=False` scores variants the GWAS never studied

**Severity: medium. Requires the user to explicitly opt out of two safe defaults**
(`split_multi=True`, `strict_allele_match=True`; `_config.py:70-72`).

The dosage arithmetic is **not** wrong — `_calculate_dosage` string-matches the
effect allele, so whatever it counts, it counts truthfully. What is lost is
**variant identity**: the guarantee that the variant at a coordinate is the variant
the weights row describes.

With `strict_allele_match=False` there is **no allele check at all** — the only
filter is `hl.is_defined(mt.weights_info)` on a locus-only join. And
`effect_allele == alleles[0]` (effect allele is REF) is satisfied by *any* variant
at that locus, because the reference base is the same string whichever ALT is there:

```
chr1:3000  A/T  -> alleles[0] = "A"
chr1:3000  A/C  -> alleles[0] = "A"
```

**Worked example.** Weights say `chr1:3000`, effect `A`, noneffect `G` — an A/G SNP.
The VDS has an **A/T** SNP at that position. The A/G SNP is not there.

| sample | genotype | copies of A counted | |
|---|---|---|---|
| S1 | A/A | — | never visited (Finding 1) |
| S2 | A/T | **1** | phantom contribution |
| S3 | T/T | 0 | |

S2 alone picks up a copy, purely because it happens to be A/T rather than T/T. The
weight was estimated for an A-vs-G contrast and is being driven by an unrelated T
genotype. Per-sample noise, so it perturbs rankings.
(`test_non_split_loose_scores_a_variant_the_gwas_never_studied`)

`strict_allele_match=True` asks whether one weights allele is REF *and* the other is
a real ALT here — `(A==A) & {T}∋G` is false, `(G==A)` is false — and drops the row
for everyone. No phantom term.

**The mirror case fails safe.** If the effect allele is an ALT that isn't present
(`chr1:4000`, effect `G`, VDS `A/T`), nothing string-matches `G`, dosage is 0, and
the row contributes nothing — though it is still counted in `n_matched`.

### Why a non-match usually means the *weights* are wrong

If a GWAS names an A/G SNP and the AoU VDS has no such variant, the likely
explanation is a **bad weights row**, not a gap in AoU — a callset that size is not
missing common GWAS SNPs. Realistic causes: **strand flip** (the library does *no*
strand harmonization anywhere; palindromic A/T and C/G variants are undetected),
**build mismatch** (GRCh37 coordinates against a GRCh38 VDS), bad rsID→allele
mapping, or a genuinely filtered variant. In every one of those cases the right
action is to drop the row — which is what strict does.

---

## Finding 4 — the split path silently drops rows of the opposite orientation

`ref_is_effect_allele` is a **global** flag, but allele orientation is a
**per-row** property of the weights file. The split path builds its `(locus,
alleles)` join key from the flag, so every weights row whose orientation disagrees
is silently dropped.

On the 5-row mock file: the flag off matches 2 rows, on matches 1. **No setting
scores them all**, and the score still comes out as a clean float.
(`test_split_silently_drops_rows_of_the_opposite_orientation`)

Whether this bites in practice depends on whether real PGS Catalog files are
uniformly oriented. Worth checking before treating it as a bug.

---

## How each config counts the effect allele

Dosage on the non-split path is always `_calculate_dosage`; on the split path it is
always `GT.n_alt_alleles()`. In **all** configs, a sample absent from `variant_data`
contributes nothing. The configs differ in which rows survive to be counted.

| config | join key | allele check | hom-ref sample | effect allele on REF |
|---|---|---|---|---|
| `split_multi=True` (default) | `(locus, alleles)` | the join key enforces it | contributes 0 (correct) | **row dropped** unless `ref_is_effect_allele=True` |
| `split_multi=True`, `ref_is_effect_allele=True` | `(locus, alleles)` | the join key enforces it | contributes 0 | scored as `-w·n_alt`; uniform `2w` offset, rankings safe |
| `split_multi=False`, `strict_allele_match=True` | locus only | `_check_allele_match` | contributes 0 (**wrong**, Finding 1) | scored, but hom-ref samples lost |
| `split_multi=False`, `strict_allele_match=False` | locus only | **none** | contributes 0 (**wrong**) | scored against whatever variant is there (Finding 3) |

`n_matched` counts rows that survived the row filter, so under the loose setting a
row that matches on locus but contributes dosage 0 still counts as "matched" — it
over-reports overlap with the VDS.

---

## Open decisions

**A. Finding 1 (hom-ref lost on the non-split REF-effect path).** Needs a real fix;
this is the only finding that silently reorders samples under a *documented*
configuration. The fix is not a one-liner: dosage must be computed for samples that
have no entry, which means either densifying against the reference blocks
(expensive) or restructuring the aggregation so absent samples get an explicit
hom-ref default. Given the split path is the default and handles this correctly,
**deprecating or restricting the non-split path is a serious alternative** to fixing
it.

**B. Finding 3 (`strict_allele_match=False`).** Three options:

1. **Targeted fix.** Require the noneffect allele to be confirmed *only when the
   effect allele is REF* — the sole case that fails dangerously — while continuing
   to tolerate an unverified noneffect allele when the effect allele is an ALT,
   which fails safe anyway. Keeps the one thing the loose setting legitimately buys
   (tolerating a junk "other allele" column) and removes the corruption.
2. **Remove `strict_allele_match` entirely**, always strict.
3. **Leave it, document it loudly.**

Option 1 remains the recommendation. Note that if A is resolved by dropping the
non-split path, B disappears with it.

**C. Finding 2 (no-call).** Confirm first whether no-call entries occur in the AoU
`variant_data` at all.

**D. Finding 4 (orientation).** Confirm first whether real weights files are
uniformly oriented.

---

## Other open items

- **`main` is behind `dev`.** All current work is on `dev`. Clean fast-forward
  whenever wanted.
- **`_standardize_chromosome_column` and scaffold contigs.** The vectorized rewrite
  (`0c6f467`) prefixes every unprefixed contig, so GRCh38 scaffold/alt contigs
  (`KI270728.1`) become `chrKI270728.1` and make `hl.locus` throw. Fails loudly, and
  PGS Catalog / PRS-CS files don't emit scaffold contigs, so this is theoretical —
  but it is a narrow regression.
- **Bundled weights data is still unused by tests.** `aoutools/data/*.{csv,tsv}`
  are 100 real GRCh38 chr1 variants whose `ID` column encodes `chr:pos:REF:ALT`.
  They would give `_reader.py` its first real-hail coverage (it is 100% mocked
  today) and exercise chunking at a more realistic scale.
