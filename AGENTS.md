# AGENTS.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

`aoutools` is a Python library for *All of Us* Researcher Workbench analysis.
The only submodule so far is `aoutools.prs`, which reads polygenic score weight
files and calculates PRS directly against the All of Us Hail VDS.

## Environment & commands

Uses **pixi** (see `pyproject.toml`) with four environments: `default`/`dev`
(development + docs), `ci` (linux-64, pins `hail==0.2.135` and the rest of the
Workbench genomics runtime), `integration` (real hail + a JDK, **both** platforms),
and `lint` (ruff only).

The `default` env is auto-activated by **direnv** (`.envrc`, committed) on `cd`
into the repo, so its tools are already on `PATH`; run `direnv allow` once per
clone. Everything else still goes through `pixi run …` — never system python.

```bash
pixi run -e ci test                    # mocked suite (tests/prs); Linux only
pixi run -e ci pytest tests/prs/test_reader.py::test_name   # single test
pixi run -e integration test-integration   # real-hail suite; macOS or Linux
pixi run docs                    # build Sphinx HTML docs
pixi run lint                    # ruff check + format --check (what CI runs)
pixi run format                  # ruff format + check --fix (writes)
pixi run setup-hooks             # install pre-commit hooks; once per clone
```

## Two test tiers

**`tests/prs/` mocks `hail`** with `MagicMock` (but still imports it, so it runs
only under `ci`/Linux, not macOS). It checks that the code calls the right hail
methods. **It cannot tell a correct score from a wrong one** — a flipped effect
allele or a bad join key still yields a clean float column and a green suite.

**`tests/integration/` runs real hail.** It starts a local Spark backend, builds a
GRCh38 mock VDS and a mock GWAS summary, and asserts **per-sample allele copy
numbers**. All weights are `1.0`, so `prs` *is* the copy count and the assertions
read as arithmetic.

The load-bearing fact it pins: **a hom-ref sample has no entry in `variant_data`**,
and hail filters absent entries out of the entry stream, so aggregators never visit
it. No missing-genotype default can reach it. Confirmed on the real All of Us VDS
(`notebooks/verify_hom_ref_dosage.ipynb`); it is why the non-split scoring path was
removed.

Several tests deliberately pin **known bugs** so a fix is reviewable; each says so
in its docstring and points at `TODO.md`. Do not "make them pass" — read them.

This tier is why `hail` is installed from **PyPI** in the `integration` feature
rather than conda: the conda package is linux-64 only, but the PyPI wheel is
`py3-none-any` (it bundles the Spark JAR), so given a JDK it runs on macOS too.
`hail` stays out of `[project.dependencies]` regardless: the Workbench provides
it, and it pins numpy/pandas tightly.

Any change to scoring logic should land with an integration test. The mocked tier
will not catch you.

## The third tier: `notebooks/`, run on the Workbench

Both offline tiers have a hole only the real Workbench can fill, so two notebooks
are checked in. **They are not run by CI** — a human runs them on a Hail Genomic
Analysis environment before a release.

- `validate_scoring_on_aou.ipynb` — the real-data counterpart of
  `tests/integration/`. Every fixture in `tests/integration/conftest.py` is a
  *claim about the shape of the real VDS*; if a claim is false the mock suite
  stays green and the scores are wrong anyway. So this notebook **finds a real
  variant of each shape** in the real VDS and checks the library against copy
  numbers counted from `hl.vds.to_dense_mt` — an oracle that never calls
  `min_rep`, `split_multi`, or `aoutools`. That independence is not ceremony:
  writing it is what turned up **Finding 6**.
- `validate_public_api_on_aou.ipynb` — `calculate_prs`, `calculate_prs_batch`,
  and `calculate_pgs` all hard-raise unless `output_path` starts with `gs://`, so
  **no offline test reaches any of them**. Nor the PGS Catalog download, nor
  `_stage_local_file_to_gcs`. This notebook is the only check on all of it.

`aoutools.init_hail()` and `aoutools.get_vds_path()` (`_workbench.py`) exist so
the Workbench boilerplate lives in one place: requester-pays billing from
`$GOOGLE_PROJECT`, the reference genome set *after* `hl.init` (passing
`default_reference` to it is deprecated), and a fallback for `$WGS_VDS_PATH`,
which current Workbench images no longer export. The fallback names a specific
All of Us data release and **warns**; bump `DEFAULT_VDS_PATH` when AoU cuts a new
one. It is deliberately not auto-discovered — picking "the newest version in the
bucket" could hand back genomic data that does not match the CDR the workspace is
registered against, which is a wrong analysis rather than an error.

Two older notebooks record findings rather than validate the library.
`verify_hom_ref_dosage.ipynb` (why the non-split path was removed) is **frozen
evidence** — its cells construct `PRSConfig(split_multi=False)` and raise
`TypeError` today, so they keep the pre-helper `hl.init` boilerplate on
purpose; modernizing cells that cannot run would only make them look as if they
could. `measure_minrep_locus_shift.ipynb` (why Finding 5 was closed) **is still
runnable**, and it is what you re-run if the locus-shift tripwire in
`validate_scoring_on_aou.ipynb` ever fires on a new VDS release.

## Lint & formatting

**ruff** is both linter and formatter; config is `[tool.ruff]` in
`pyproject.toml` (`E,W,F,I,UP,B` at 80 cols). Pre-commit runs both on commit.

Sphinx **mocks `hail`** when building docs, and a mocked `hl.Table` does not
support the PEP 604 `|` operator. Modules that annotate with hail types in a
union (`_config.py`, `_calculator_utils.py`) therefore need
`from __future__ import annotations` so those annotations are never evaluated at
import — without it, autodoc fails to import the module and **silently drops the
whole API reference** while still exiting 0. Keep that import when adding
hail-typed unions elsewhere.

## Architecture

Public API is re-exported from internal `_`-prefixed modules via
`aoutools/prs/__init__.py`:

- `read_prs_weights`, `read_prscs` (`_reader.py`) — parse weight files (headers
  optional, custom `column_map`, allele validation) into Hail Tables.
- `calculate_prs` (`_calculator.py`) — score one weights table against a VDS.
- `calculate_prs_batch` (`_calculator_batch.py`) — score many at once, filtering
  the VDS once against the union of all loci.
- `calculate_pgs` (`_workflow.py`) — end-to-end: download PGS Catalog files
  (explicit PGS IDs only, not EFO/PGP, to avoid unbounded downloads), read, batch.
- `download_pgs` (`_downloader.py`), `PRSConfig` (`_config.py`, all calc params).

**Calculation strategy** (the core design): to stay cheap on the large VDS, the
weights table is chunked into `chunk_size` variants; each chunk builds 1bp
intervals for `hl.vds.filter_intervals` so only relevant loci are read, computes
its PRS, and results are summed across chunks. There is **one** scoring path:
`hl.vds.split_multi` splits multi-allelic sites, then the weights are matched to
each row. Shared helpers live in `_calculator_utils.py`.

**The join is by locus, then alleles locally — never a re-key.** The split
`variant_data` is keyed by `(locus, alleles)` in the VDS's own — possibly
non-minimal — representation, while a GWAS names variants minimally. Matching on
the full row key would therefore mean rewriting the VDS alleles to minimal form,
i.e. `key_rows_by`, which **shuffles the whole chunk's entries**. Instead
`_group_weights_by_locus` groups the (small) weights into one array per locus;
`weights_by_locus[mt.locus]` is a key-*prefix* join against the MT's leading
key, which hail does with no shuffle; and `_match_weight_at_locus` then picks the
row whose **unordered allele set** equals the row's minimal `canonical_alleles`
(annotated by `_split_multi_with_total_dosage`). Only the small weights table
moves. The set match is orientation-agnostic by construction, so it also handles
a file that writes the pair either way round.

Two silent-drop bugs are pinned here, both invisible until real data. A `G/A`
SNP where REF sorts after ALT (`chr1:8000`) broke an earlier `hl.sorted` key
that only matched REF-before-ALT variants — it cost ~half of every real score
(922 of 1,940 for PGS000746). A non-minimal biallelic like `AAAG/GAAG`
(`chr1:8500`), which `split_multi` passes through un-normalized, never matched
its minimally-named `A/G` weights. The allele-set-on-`canonical_alleles` match
closes both. Do **not** "simplify" this to a keyed `(locus, alleles)` join: it
reintroduces one or both drops, or the entry-shuffling re-key.

**The hom-ref offset.** Orientation is resolved per row against the VDS's REF.
When the effect allele is the REF, the true contribution is `w·(2 − n_non_ref)`,
which splits into a per-entry `−w·n_non_ref` and a row-level `2w`. The entry term
is already right for hom-ref samples (`n_non_ref` is 0); the `2w` is added as a
**row-level constant**, because a hom-ref sample has no entry and no entry
aggregator can reach it. Three invariants hold it up, all of them silent if
broken: the offset is only summed over rows that *matched* a variant, it is
cancelled for no-call entries (`_entry_contribution`), and it never varies per
sample. `tests/integration/` asserts each.

**Two dosages, because splitting downcodes.** `split_multi` rewrites every *other*
ALT to the reference, so at the `[C,T]` row of a `C/G,T` site a C/G carrier has
`GT == 0/0` — indistinguishable from a hom-ref sample. That is correct for an
ALT-effect weight (it carries no T) and **wrong** for a REF-effect one (it carries
one C, not two). So `_split_multi_with_total_dosage` annotates each entry with
`n_non_ref` — the total non-ref count from the **pre-split** local genotype, which
survives the split — and `_entry_contribution` counts `n_non_ref` for REF-effect
rows and `GT.n_alt_alleles()` for ALT-effect rows. Both branches are load-bearing;
using either one everywhere turns integration tests red. This was Finding 6.

Two knobs were removed for silently losing data: `split_multi=False` selected a
path that zeroed every hom-ref sample at a REF-effect variant (reordering the
cohort), and `ref_is_effect_allele` declared orientation file-wide when it is a
per-row property, dropping every row that disagreed. `TODO.md` has the evidence.

**`hl.vds.split_multi` is called with the default `filter_changed_loci=False`,
which raises. That is deliberate — don't set it to `True` to silence an error.**
`min_rep` trimming a shared *suffix* is safe and relied on (`[AGGGC, A, GGGGC]` →
`A/G` at the same locus, which is how a GWAS names it). Trimming a shared *prefix*
would *move* the locus, and hail can then only raise or silently drop the allele.
No such variant exists in AoU (0 of 6,001,424 ALTs measured), so the exception is
a tripwire for a future VDS release, not a crash risk.

**`split_multi` only min_reps the rows it actually splits**, which is why
`_split_multi_with_total_dosage` annotates a `canonical_alleles` row field (the
`min_rep` of every row) for the join above to match on. An already-biallelic row
is passed through with its *original* alleles, so without this a non-minimal
biallelic variant — `AAAG/GAAG`, an A→G SNP a GWAS names `A/G` — keeps
`AAAG/GAAG` and never joins the minimally-named weights, scoring 0 with no error.
`min_rep` is idempotent on the rows split_multi already reduced, so annotating
every row only changes the passthroughs, and it is annotate-not-re-key precisely
so the entry data never shuffles. The locus-shift tripwire (`or_error` if
`min_rep` moves the locus) is carried onto these passthrough rows too, which
`filter_changed_loci` does not cover. The multi-allelic non-minimal case is
normalized as a *side effect* of splitting and so hid this for a long time;
`validate_scoring_on_aou.ipynb` found it on the real VDS (`chr1:1409159
AAAG/GAAG`) because its non-minimal fixture is biallelic where
`tests/integration/`'s was only multi-allelic. Both shapes are now pinned.

`_utils.py:_stage_local_file_to_gcs` copies local paths into
`$WORKSPACE_BUCKET/data/...` because Hail's Spark cluster can't read the local
notebook filesystem. Package version derives from installed metadata
(`pyproject.toml` is the single source of truth); don't hardcode it.
