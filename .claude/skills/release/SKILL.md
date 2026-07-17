---
name: release
description: Cut and publish an aoutools release to PyPI — bump the version, build, upload with twine, and tag. User-invoked only; a PyPI upload cannot be undone.
disable-model-invocation: true
---

# Release aoutools to PyPI

There is no publish workflow — releases are built and uploaded by hand from a
local checkout. **A version uploaded to PyPI can never be replaced or re-used**,
even after deletion. Get it right before uploading, not after.

Ask the user for the new version if they did not say. Current released versions:
0.1.0, 0.1.1, 0.1.2.

## 1. Preconditions

- Working tree clean, on `main`, with `dev` already merged in.
- CI green on that commit (`gh run list --branch main --limit 1`).
- Run `/verify` — lint, docs, and the test suite must all pass.

## 2. Bump the version

`pyproject.toml` is the **only** place the version appears:

```toml
version = "X.Y.Z"
```

Do not add it anywhere else. `aoutools/__init__.py` and `docs/source/conf.py`
both read it from installed package metadata (`importlib.metadata`) precisely so
it cannot drift — `conf.py` used to hardcode it and silently published the wrong
version for two releases.

Then re-sync the editable install so the metadata reflects the bump:

```bash
pixi install
pixi run -e default python -c "import aoutools; print(aoutools.__version__)"
```

Commit the bump.

## 3. Clean dist/ — do not skip this

```bash
rm -rf dist/
```

`dist/` is gitignored and is **not** cleaned by the build. It currently holds
artifacts from earlier releases. If you build without clearing it, `twine upload
dist/*` picks up the old files too, PyPI rejects the already-published version,
and the whole upload fails partway. Always start from an empty `dist/`.

## 4. Build and check

```bash
pixi run -e default python -m build     # python-build is in the dev feature
pixi run -e default twine check dist/*
```

Confirm `dist/` contains exactly one wheel and one sdist, both at the new
version, and nothing else.

## 5. Upload

Irreversible. Confirm with the user before running:

```bash
pixi run -e default twine upload dist/*
```

## 6. Tag

The scheme is a `v` prefix: `v0.1`, `v0.1.1`.

```bash
git tag vX.Y.Z
git push origin vX.Y.Z
```

**Note:** `v0.1.2` was published to PyPI on 2026-02-06 but originally shipped
untagged. The tag has since been backfilled locally at `b03c522` (the commit that
bumped the version); confirm it is pushed to the remote (`git push origin
v0.1.2`) if the release predates you.

## 7. Confirm

- PyPI shows the new version: `https://pypi.org/project/aoutools/`
- Read the Docs built the new tag and shows the right version in the sidebar —
  it reads the version from package metadata, so it should now match `pyproject.toml`.
