# AGENTS.md

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
(`notebooks/validate_scoring_on_aou.ipynb`, check 1); it is why the non-split
scoring path was removed.

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

Offline tiers have holes only the real Workbench can fill. **Not run by CI** — a
human runs these on a Hail Genomic Analysis environment before a release.

- `validate_scoring_on_aou.ipynb` — the real-data counterpart of
  `tests/integration/`. Each `conftest.py` fixture is a *claim about the shape of
  the real VDS*; if a claim is false the mock suite stays green and the scores are
  wrong anyway. So this **finds a real variant of each shape** and checks the
  library against copy numbers from `hl.vds.to_dense_mt` — an oracle that never
  calls `min_rep`, `split_multi`, or `aoutools`. Writing it turned up **Finding
  6**. Check 1 pins the hom-ref-has-no-entry fact on real data.
- `validate_public_api_on_aou.ipynb` — `calculate_prs`, `calculate_prs_batch`,
  and `calculate_pgs` all hard-raise unless `output_path` starts with `gs://`, so
  **no offline test reaches them**, nor the PGS Catalog download nor
  `_stage_local_file_to_gcs`. The only check on all of it.
- `validate_synthetic_control_on_aou.ipynb` — the **positive-control** tier.
  Rather than *finding* a real variant of each shape, it **builds** a synthetic
  VDS and weights file in which every scoring path is present by construction,
  computes the expected PRS independently (pure Python over the design,
  cross-checked against a `to_dense_mt` oracle), and drives the **public API**
  against that control — so it checks the user-facing functions and the `gs://`
  round-trip on a *known answer*, which `validate_public_api_on_aou.ipynb` (real
  PGS files, no ground truth) cannot. Ships a diagnostic toolkit that decomposes
  a discrepancy per variant and classifies its signature (constant offset /
  hom-ref-only / genotype-dependent / whole-variant drop) onto the matching
  `TODO.md` finding. Because its data is synthetic real-hail (not the real
  VDS), every cell except the `gs://` public-API ones runs offline under
  `pixi run -e integration` too — the same tier `tests/integration/` uses.
- `measure_minrep_locus_shift.ipynb` — records why Finding 5 is closed (locus-
  shift rate is zero in AoU); re-run it if the tripwire below ever fires.

`aoutools.init_hail()` and `aoutools.get_vds_path()` (`_workbench.py`) keep the
Workbench boilerplate in one place: requester-pays from `$GOOGLE_PROJECT`, the
reference set *after* `hl.init` (passing `default_reference` to it is deprecated),
and a fallback for `$WGS_VDS_PATH`, which current images no longer export. The
fallback names a specific AoU release and **warns**; bump `DEFAULT_VDS_PATH` when
AoU cuts a new one. It is deliberately not auto-discovered — "newest in the
bucket" could return data not matching the workspace's CDR, a wrong analysis
rather than an error.

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

## Logging

Every module logs through `logging.getLogger(__name__)`, so all output lives
under the `aoutools.*` tree; `aoutools/__init__.py` attaches a `NullHandler` so
the library is silent until the application configures a handler. A notebook
opts in with `logging.basicConfig(level=...)` plus
`logging.getLogger("aoutools").setLevel(logging.INFO)`.

**Levels have fixed meanings** — hold the line on them:

- **DEBUG** — fine-grained internal steps and subprocess I/O (which validator is
  running, a command's stdout). Off by default; for diagnosing, not narrating.
- **INFO** — user-meaningful milestones only: a file loaded and its variant
  count, per-chunk progress, a stage completing with its timing. Default-INFO
  output should read as a short progress story, not a trace.
- **WARNING** — data-quality events the user should notice but that do not stop
  the run: variants dropped for bad alleles or missing coordinates, a scoring
  file skipped, an empty result. Routine, by-design drops (zero-weight rows) are
  INFO, not WARNING — do not cry wolf.
- **ERROR** — reserved for a failure that is **not** also raised. Never
  `logger.error(...)` right before `raise`: the exception already carries it,
  and with the `NullHandler` the log line is invisible anyway. Put failure
  detail in the exception message.

**Two channels, one rule.** `warnings.warn` is for actionable *setup* advisories
aimed at the user (a fallback VDS path, no billing project, a recovered
bucket) — things they should change, shown once. `logger.*` is for runtime
narration of a call that is proceeding normally. `_workbench.py` uses
`warnings.warn`; everywhere else uses the logger.

**Style.** Noun-phrase milestone or past-tense completion; no trailing `...`
(the one exception is `_log_timing`, whose start line is a deliberate
in-progress marker). Always use lazy `%s`/`%d` args, never an f-string inside
the log call. Keep the wording plain — the audience is researchers, not
engineers (see the `docs/source/` voice).

There is **no per-call verbosity flag.** A user raises detail by lowering the
`aoutools` logger level, not by passing `verbose=`. `PRSConfig.detailed_timings`
is the one exception and is scoped to timing granularity, not a general switch.

## Architecture

Public API is re-exported from internal `_`-prefixed modules via
`aoutools/prs/__init__.py`:

- `read_prs_weights` (`_reader.py`) — parse weight files (headers optional,
  custom `column_map`, allele validation) into Hail Tables. `read_prscs` is a
  **deprecated** thin wrapper (fixed PRS-CS column layout); it warns and defers
  to `read_prs_weights` — don't build on it.
- `calculate_prs` (`_calculator.py`) — score one weights table against a VDS.
- `calculate_prs_batch` (`_calculator_batch.py`) — score many at once, filtering
  the VDS once against the union of all loci.
- `calculate_pgs` (`_workflow.py`) — end-to-end: download PGS Catalog files
  (explicit PGS IDs only, not EFO/PGP, to avoid unbounded downloads), read, batch.
- `download_pgs` (`_downloader.py`), `PRSConfig` (`_config.py`, all calc params).

**Calculation strategy** (the core design): to stay cheap on the large VDS, the
weights are chunked into `chunk_size` variants; each chunk builds 1bp intervals
for `hl.vds.filter_intervals` so only relevant loci are read, and results are
summed across chunks. One scoring path: `hl.vds.split_multi` splits multi-allelic
sites, then the weights are matched to each row. Helpers live in
`_calculator_utils.py`.

**The join must survive two mismatches between the weights and the VDS, each of
which was a silent-drop bug.** (1) *Allele order* — orientation is per-row, not
file-wide (the PGS Catalog does not harmonize the effect allele onto the ALT), so
a `G/A` SNP whose REF sorts after its ALT (`chr1:8000`) must still match. An
earlier `hl.sorted` key matched only REF-before-ALT variants and dropped ~half of
every score (922 of 1,940 for PGS000746). (2) *Minimal representation* —
`split_multi` only `min_rep`s the rows it splits, so an already-biallelic
non-minimal variant (`AAAG/GAAG`, a GWAS's `A/G`, `chr1:8500`) passes through
un-normalized and never matches. Both are closed together:
`_split_multi_with_total_dosage` annotates `canonical_alleles` (the `min_rep` of
every row; idempotent on split rows), and `_match_weight_at_locus` matches on the
unordered allele **set** against it. Both shapes are pinned in `tests/integration/`.

The match is also **shuffle-free**, which is why it is not a keyed join.
`_group_weights_by_locus` groups the small weights into one array per locus;
`weights_by_locus[mt.locus]` is a key-*prefix* join (no shuffle), then alleles
match locally. Keying on `(locus, alleles)` instead would force a `key_rows_by`
to minimal alleles, **shuffling every entry**. Do **not** "simplify" to a keyed
`(locus, alleles)` join: it reintroduces a silent drop or the re-key shuffle.

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

**The locus-shift tripwire.** `min_rep` trimming a shared *suffix* is safe and
relied on (`[AGGGC, A, GGGGC]` → `A/G`, same locus, how a GWAS names it). Trimming
a shared *prefix* would *move* the locus, which hail can only raise on or silently
drop. `hl.vds.split_multi` runs with the default `filter_changed_loci=False`
(raises), and the `canonical_alleles` annotation carries the same `or_error` guard
onto biallelic passthroughs that split_multi doesn't cover. No such variant exists
in AoU (0 of 6,001,424 ALTs), so this is a tripwire for a future VDS release, not
a crash risk. The raise is a **split-step** guard, so it only catches a shifted
variant that reaches `split_multi`: the per-chunk interval prefilter is built at
the *weights* locus, so `filter_intervals` drops a variant whose row sits upstream
of its minrep'd, GWAS-named locus *before* it is split — making a downstream-named
shift a silent `n_matched` shortfall, not a raise. `tests/integration/` and
`notebooks/validate_synthetic_control_on_aou.ipynb` pin the raise through the
`_calculate_prs_chunk` seam (which splits the unfiltered VDS) for that reason.
Don't set `filter_changed_loci=True` to silence it.

`_utils.py:_stage_local_file_to_gcs` copies local paths into
`$WORKSPACE_BUCKET/data/...` because Hail's Spark cluster can't read the local
notebook filesystem. Package version derives from installed metadata
(`pyproject.toml` is the single source of truth); don't hardcode it.
