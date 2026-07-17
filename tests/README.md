# Tests

Two tiers, run separately:

```bash
pixi run -e ci test                    # tests/prs/ -- mocks hail; Linux only
pixi run -e integration test-integration   # tests/integration/ -- real hail; macOS or Linux
```

- **`tests/prs/`** mocks `hail`, so it checks that the code *calls* the right hail
  methods — it **cannot tell a correct score from a wrong one**. A change to
  scoring logic that passes here can still be wrong.
- **`tests/integration/`** runs real hail against a small mock VDS and asserts
  per-sample allele copy numbers. This is the tier that catches a wrong score, so
  any scoring change should land with a test here.

> [!WARNING]
> Some tests **deliberately pin a known bug** so a future fix is reviewable. Each
> says so in its docstring and points at `TODO.md`. If one of these fails, do
> **not** "make it pass" — read the docstring first.

For why the tiers are split this way (and why `hail` comes from PyPI, not conda),
see the **"Two test tiers"** section of `AGENTS.md`.
