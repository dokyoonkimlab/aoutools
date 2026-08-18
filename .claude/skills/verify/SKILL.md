---
name: verify
description: Verify a change to aoutools before committing — lint, docs, and the full test suite, including how to run the hail-dependent tests on macOS. Use after changing any code in this repo, or when asked to check that a change works.
---

# Verify an aoutools change

Run these in order. They are fast — the whole set is under a minute.

## 1. Lint and format (this is exactly what CI runs)

```bash
pixi run lint      # ruff check + ruff format --check; never writes
pixi run format    # ruff format + check --fix; writes
```

## 2. Docs

```bash
pixi run docs
```

This builds with `-W` (warnings-as-errors, set in `docs/Makefile`). **Do not
relax that.** Sphinx reports a failed autodoc import as a *warning* and then
renders an empty page, so without `-W` the build exits 0 while silently
publishing an API reference with no content. That has already happened once.

Sphinx **mocks `hail`**, and a mocked `hl.Table` does not support the PEP 604
`|` operator. Any module that annotates with a hail type in a union needs
`from __future__ import annotations` (see `_config.py`, `_calculator_utils.py`),
or the module fails to import under autodoc and the API reference vanishes.

## 3. Tests

`hail` has no macOS wheel, so `pixi run -e ci test` cannot run on a Mac —
that environment is linux-64 only. **But the mocked suite only needs hail to be
importable**, so a stub is enough to run that tier locally:

```bash
STUB=$(mktemp -d)
mkdir -p "$STUB/hail/vds" "$STUB/hail/utils" "$STUB/hailtop"
printf 'from unittest.mock import MagicMock\nclass Table: pass\nclass MatrixTable: pass\ndef __getattr__(n): return MagicMock()\n' > "$STUB/hail/__init__.py"
printf 'from unittest.mock import MagicMock\nclass VariantDataset: pass\ndef __getattr__(n): return MagicMock()\n' > "$STUB/hail/vds/__init__.py"
printf 'from unittest.mock import MagicMock\ndef __getattr__(n): return MagicMock()\n' > "$STUB/hail/utils/__init__.py"
printf 'class FatalError(Exception): pass\nclass HailUserError(Exception): pass\n' > "$STUB/hail/utils/java.py"
printf '' > "$STUB/hailtop/__init__.py"
printf 'from unittest.mock import MagicMock\ndef __getattr__(n): return MagicMock()\n' > "$STUB/hailtop/fs.py"

PYTHONPATH="$STUB" pixi run -e default python -m pytest tests -q -m 'not integration'
```

`-m 'not integration'` is **required**, and is what makes this match CI's
`pixi run -e ci test`. Without it pytest also collects `tests/integration/`,
which drives a real Spark backend and cannot run against a stub — you get errors
that look like a broken suite but only mean the wrong tier was selected.

The stub must cover every hail name the package imports at module scope, not just
the ones the tests touch. `hail.utils.java.FatalError` is one such import
(`_calculator_utils.py`), and a missing name fails at **collection**, before any
test runs. If you see `ModuleNotFoundError: No module named 'hail.<something>'`,
add it here rather than assuming the suite is broken.

Expect **78 passed, 26 deselected** in about a second, matching CI.

**The stub is a convenience, not the source of truth.** It is not hail: it cannot
catch anything that depends on real hail behavior — typechecks on expressions,
`hl` constructors rejecting a bad argument, or an import that only breaks against
the real package. A green stub run is necessary, not sufficient. On Linux, or in
CI, run the real thing:

```bash
pixi run -e ci test        # Linux only; real hail 0.2.135
```

CI runs three jobs (`lint`, `docs`, `test`) on every push to `main`/`dev`. It is
the authority.

## 4. If you touched the scoring path

Changes to `_calculator.py`, `_calculator_batch.py`, `_calculator_utils.py`,
`_reader.py`, or `_config.py` can produce a **wrong-but-passing** score: the
tests mock hail, so a flipped effect allele or a bad join key still yields a
clean float column and a green suite. Run the `prs-correctness-reviewer`
subagent on the diff. Do not treat a passing test run as evidence the score is
right.
