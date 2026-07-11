# TODO

## Open decision: what to do about `strict_allele_match=False`

**Status:** undecided. No code has been changed. `aoutools/prs/` is untouched —
so far this is analysis only, plus a correction to the `prs-correctness-reviewer`
brief, which had described this behavior wrongly.

**Reachability:** the defaults are safe. `split_multi=True` and
`strict_allele_match=True` (`_config.py:70-72`). A user has to explicitly set
*both* off to reach the behavior below.

### How each config counts the effect allele

Dosage is always computed by `_calculate_dosage`
(`_calculator_utils.py:214-280`), which **string-matches the effect allele**
against the sample's reconstructed genotype and counts copies. It never indexes a
position blindly. The three configs differ only in **which rows survive to reach
that counting step**.

| config | join key | allele check before counting |
|---|---|---|
| `split_multi=True` (default) | `(locus, alleles)` | none needed — the join key itself enforces it |
| `split_multi=False`, `strict_allele_match=True` (default when non-split) | locus only | `_check_allele_match` drops rows whose alleles don't correspond to the variant at that position |
| `split_multi=False`, `strict_allele_match=False` | locus only | **none** |

Missing genotypes are treated as homozygous reference (dosage `2` if the effect
allele is REF, else `0`). This is only correct because the VDS stores hom-ref
sparsely: at a given row, `variant_data` contains entries **only** for samples
carrying a non-reference call. Everyone else is absent, not `GT=[0,0]`.

### The problem, precisely

The dosage arithmetic is **not** wrong. Whatever it counts, it counts truthfully.
What `strict_allele_match=False` loses is **variant identity** — the guarantee
that the variant sitting at a coordinate is the variant the weights row describes.

The trap is that `effect_allele == alleles[0]` (effect allele is the REF base) is
satisfied by *any* variant at that locus, because the reference base is the same
string no matter which ALT is present:

```
chr1:1000  A/T   → alleles[0] = "A"
chr1:1000  A/C   → alleles[0] = "A"
chr1:1000  A/AT  → alleles[0] = "A"
```

So a REF-effect-allele weights row matches on the coordinate alone. The check
carries no information about whether the right variant was found.

### Worked example

GWAS row: `chr1:1000`, effect allele `A`, noneffect allele `G` (so the GWAS
studied an **A/G** SNP, and `A` is also the reference base).

VDS at `chr1:1000`: the A/G SNP is **not there**. A different SNP occupies the
position — **A/T**, so `alleles = ["A", "T"]`.

The locus-only join attaches the A/G weights row to the A/T variant. Nothing
compares `G` to `T`.

**`strict_allele_match=False`** — no check runs; the row is counted:

| sample | bases at chr1:1000 | in `variant_data`? | copies of `A` counted |
|---|---|---|---|
| S1 | A/A | absent (sparse hom-ref) | 2 |
| S2 | A/T | `GT=[0,1]` | 1 |
| S3 | T/T | `GT=[1,1]` | 0 |

Every count is a *truthful* count of `A` — S2 really does carry one `A`. The
defect is that the score now contains a term that varies with each sample's **T**
genotype, an allele the GWAS never studied, weighted by an effect size estimated
for an A-vs-G contrast. Two people who are identical with respect to the GWAS
variant get different scores. The effect is **per-sample**, so it perturbs
rankings, percentiles and standardized scores — not just absolute values.

**`strict_allele_match=True`** — `_check_allele_match` asks whether one weights
allele is REF *and* the other is a real ALT here:

```
("A" == REF "A")  &  {"T"} contains "G"   → True  & False → False
("G" == REF "A")  &  {"T"} contains "A"   → False & False → False
                                          → row dropped
```

The variant is excluded for everyone. The score simply loses that SNP — no
phantom term, and every sample is affected identically.

**Mirror case fails safe.** If the effect allele had been an ALT (`G`) rather than
REF, the string match against `["A","T"]` finds nothing, dosage is `0`, and the
row contributes nothing to anyone even under the loose setting.

### Why a non-match usually means the *weights* are wrong

If a GWAS names an A/G SNP at a position and the All of Us VDS has no such
variant, the likely explanation is a **bad weights row**, not a gap in AoU — a
callset of that size is not missing common GWAS SNPs. Realistic causes:

- **strand flip** (the library does *no* strand harmonization anywhere;
  palindromic A/T and C/G variants are not detected)
- **build mismatch** (GRCh37 coordinates against a GRCh38 VDS — the position then
  points at an unrelated variant)
- bad rsID→allele mapping in whatever produced the weights file
- genuinely rare / uncalled variant, filtered from the callset

In every one of those cases the correct action is to **drop the row**. That is
exactly what `strict_allele_match=True` does. The loose setting scores it anyway
— and the more likely a row is erroneous, the more confidently it gets scored.

### What the loose setting legitimately buys

One thing, and it is narrow: it tolerates an unreliable **noneffect** column. If
the effect allele is a genuine ALT in the VDS but the weights file's "other
allele" is wrong or junk, strict drops the variant while loose counts the effect
allele correctly. Common enough in hand-built weight files.

### Options (pick one)

1. **Targeted fix — keep the escape hatch, close the hole.** Require the noneffect
   allele to be confirmed *only when the effect allele is REF* (the sole case that
   fails dangerously); keep tolerating an unverified noneffect allele when the
   effect allele is an ALT (which fails safe anyway). Preserves the legitimate use
   above, removes the corruption.
2. **Remove `strict_allele_match` entirely**, always doing the strict check. The
   loose mode's only benefit is tolerating a bad noneffect column, which is
   arguably a weights-file problem the user should fix upstream.
3. **Leave the behavior, document it loudly** — release note + docstring warning
   that `strict_allele_match=False` can inject genotype-dependent noise from
   unrelated variants.

Option 1 is the recommendation. Whichever is chosen, note that **the non-split
path currently has no test coverage for allele matching at all** — any change
here should land with tests on `_check_allele_match` and on
`_prepare_mt_non_split`'s filter behavior.

### Also worth deciding separately

`n_matched` counts rows that survived the filter. Under the loose setting, rows
that match on locus but contribute dosage `0` (effect allele absent) are still
counted as "matched", so `n_matched` over-reports overlap with the VDS.

## Other open items

- **`main` is behind `dev`.** All current work — reviewer subagent, ruff pass,
  `/verify` and `/release` skills, direnv, docs fixes — is on `dev` only. Clean
  fast-forward whenever wanted.
- **No-call genotypes are scored as hom-ref.** A genuinely no-called (`./.`)
  genotype in `variant_data` lands in the same missing branch as a sparse hom-ref
  and is scored as two copies of the effect allele when that allele is REF, rather
  than being excluded. Pre-existing, independent of `strict_allele_match`, and
  applies on the **default** config too. Worth confirming against a real VDS
  whether such entries actually occur.
- **`_standardize_chromosome_column` and scaffold contigs.** The vectorized
  rewrite (`0c6f467`) prefixes every unprefixed contig, so GRCh38 scaffold/alt
  contigs (`KI270728.1`, `GL000009.2`) become `chrKI270728.1` and make `hl.locus`
  throw. Fails loudly, and PGS Catalog / PRS-CS files don't emit scaffold contigs,
  so this is theoretical — but it is a narrow regression from the previous
  behavior.
