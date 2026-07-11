# AGENTS.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

`aoutools` is a Python library for *All of Us* Researcher Workbench analysis.
The only submodule so far is `aoutools.prs`, which reads polygenic score weight
files and calculates PRS directly against the All of Us Hail VDS.

## Environment & commands

Uses **pixi** (see `pyproject.toml`) with three environments: `default`/`dev`
(development + docs), `ci` (linux-64, pins `hail==0.2.135` and the rest of the
Workbench genomics runtime), and `lint` (ruff only).

The `default` env is auto-activated by **direnv** (`.envrc`, committed) on `cd`
into the repo, so its tools are already on `PATH`; run `direnv allow` once per
clone. Everything else still goes through `pixi run …` — never system python.

```bash
pixi run -e ci test              # full test suite (pytest -v); Linux only
pixi run -e ci pytest tests/prs/test_reader.py::test_name   # single test
pixi run docs                    # build Sphinx HTML docs
pixi run lint                    # ruff check + format --check (what CI runs)
pixi run format                  # ruff format + check --fix (writes)
pixi run setup-hooks             # install pre-commit hooks; once per clone
```

Tests **mock `hail`** but still import it, so they run only under `ci` (Linux),
not macOS; CI runs the suite on `ubuntu-latest` for `main`/`dev`. `hail` is
deliberately kept out of `[project.dependencies]`: the Workbench provides it, it
ships no wheel on some platforms (e.g. `osx-arm64`), and it pins numpy/pandas
tightly.

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
its PRS, and results are summed across chunks. `PRSConfig.split_multi` selects
two allele-handling paths — `True` (default) splits multi-allelic sites and
joins on (locus, alleles), with `ref_is_effect_allele` orienting the effect
allele; `False` keeps them, with `strict_allele_match` controlling match rigor.
Shared helpers live in `_calculator_utils.py`.

`_utils.py:_stage_local_file_to_gcs` copies local paths into
`$WORKSPACE_BUCKET/data/...` because Hail's Spark cluster can't read the local
notebook filesystem. Package version derives from installed metadata
(`pyproject.toml` is the single source of truth); don't hardcode it.
