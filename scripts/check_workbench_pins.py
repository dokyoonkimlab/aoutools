#!/usr/bin/env python
"""
Compare the Workbench pins in ``pyproject.toml`` against a real ``pip freeze``.

The ``ci`` and ``integration`` pixi environments exist to be a *replica* of the
All of Us Workbench genomics runtime: green CI is only meaningful if the
versions underneath it are the versions researchers actually run. Nothing
enforces that, because CI cannot reach a Workbench VM. Somebody has to open a
Hail Genomic Analysis environment, run ``pip freeze``, and check.

This script is the tedious half of that check. Give it the freeze file and it
reports which pins still match, which have drifted, and what the Workbench is
carrying for a few packages deliberately left unpinned.

Usage
-----
On the Workbench, in a notebook or terminal::

    pip freeze > py_pip_freeze.txt

Then locally::

    pixi run check-pins path/to/py_pip_freeze.txt
    pixi run check-pins path/to/py_pip_freeze.txt --python-version 3.11.8

``pip freeze`` does not record the Python version, so pass ``--python-version``
(read off ``python --version`` on the same machine) to check ``requires-python``
too. Without it that one check is skipped.

Exits non-zero if any pin has drifted, so it can gate a release checklist.
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path

from packaging.requirements import Requirement
from packaging.version import InvalidVersion

# Tables in pyproject.toml whose versions claim to mirror the Workbench. Each
# is checked for an exact match. Other pixi tables (dev tooling, the lint
# feature) describe our own environment and are none of this script's business.
MIRROR_TABLES = [
    (
        "tool",
        "pixi",
        "feature",
        "ci",
        "target",
        "linux-64",
        "pypi-dependencies",
    ),
    ("tool", "pixi", "feature", "ci", "target", "linux-64", "dependencies"),
    ("tool", "pixi", "feature", "integration", "pypi-dependencies"),
]

# Entries inside those tables that are not Python packages, or are our own test
# tooling rather than part of the runtime being mirrored.
NOT_FROM_THE_WORKBENCH = {"openjdk", "pytest", "pytest-mock"}

# Packages we knowingly do not pin, listed so a freeze check still surfaces
# them. protobuf is the notable one: the Workbench carries a version that
# google-cloud-storage's own metadata forbids, so no solver can reproduce it
# (see the comment in pyproject.toml). Purely informational -- never drift.
WATCHLIST = [
    "protobuf",
    "typer",
    "aiohttp",
    "orjson",
    "requests",
    "bokeh",
    "plotly",
    "rich",
]


def normalize(name: str) -> str:
    """Fold a package name to its comparable form (PEP 503)."""
    return re.sub(r"[-_.]+", "-", name).lower()


def parse_freeze(text: str) -> dict[str, str | None]:
    """
    Read ``pip freeze`` output into ``{normalized name: version or None}``.

    Real freezes from the Workbench are messier than ``name==version``. Three
    other shapes appear and all carry a version worth recovering:

    - a local wheel, ``hail @ file:///.../hail-0.2.135-py3-none-any.whl#sha...``
    - a conda-built package, ``affine @ file:///.../work`` -- no version in the
      path at all, recorded as None so it reads as "unknown", not "missing"
    - an editable install, whose version is only in the comment pip writes
      above it: ``# Editable install with no version control (pyspark==3.5.3)``
    """
    versions: dict[str, str | None] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        editable = re.match(
            r"#\s*Editable install.*\(([A-Za-z0-9_.-]+)==([^)]+)\)", line
        )
        if editable:
            versions[normalize(editable.group(1))] = editable.group(2)
            continue
        if line.startswith("#") or line.startswith("-"):
            continue

        if "==" in line:
            name, _, version = line.partition("==")
            versions[normalize(name)] = version.strip()
            continue

        if " @ " in line:
            name, _, url = line.partition(" @ ")
            versions[normalize(name)] = _version_from_url(url.strip())
    return versions


def _version_from_url(url: str) -> str | None:
    """Recover a version from a direct-reference URL, if it names one."""
    wheel = re.search(r"/([A-Za-z0-9_.-]+)-([0-9][^-]*)-py[^/]*\.whl", url)
    return wheel.group(2) if wheel else None


def pinned_version(spec: str | dict) -> str | None:
    """
    Extract the exact version from a pixi dependency spec.

    Returns None for anything that is not an exact ``==`` pin -- a range like
    ``11.*`` is not claiming to mirror anything, so there is nothing to check.
    """
    if isinstance(spec, dict):
        spec = spec.get("version", "")
    spec = str(spec).strip()
    return spec[2:].strip() if spec.startswith("==") else None


def collect_mirror_pins(config: dict) -> dict[str, str]:
    """Gather every exact pin from the tables that mirror the Workbench."""
    pins: dict[str, str] = {}
    for path in MIRROR_TABLES:
        table = config
        for key in path:
            table = table.get(key, {}) if isinstance(table, dict) else {}
        for name, spec in table.items():
            if normalize(name) in NOT_FROM_THE_WORKBENCH:
                continue
            version = pinned_version(spec)
            if version is not None:
                pins[normalize(name)] = version
    return pins


def _row(name: str, expected: str, found: str, status: str) -> str:
    return f"  {name:<24} {expected:<12} {found:<12} {status}"


def check_mirror_pins(
    pins: dict[str, str], freeze: dict[str, str | None]
) -> int:
    """Report each exact pin against the freeze. Returns the drift count."""
    print("\nPinned mirrors (must match the Workbench exactly)")
    print(_row("package", "pinned", "workbench", ""))
    drifted = 0
    for name in sorted(pins):
        expected = pins[name]
        if name not in freeze:
            print(_row(name, expected, "-", "ABSENT from freeze"))
            drifted += 1
        elif freeze[name] is None:
            print(_row(name, expected, "?", "version not recorded"))
        elif freeze[name] == expected:
            print(_row(name, expected, freeze[name], "ok"))
        else:
            print(_row(name, expected, freeze[name] or "?", "DRIFTED"))
            drifted += 1
    return drifted


def check_runtime_floors(config: dict, freeze: dict[str, str | None]) -> int:
    """
    Report the ``[project.dependencies]`` floors against the freeze.

    These are lower bounds, not pins, so the only failure is a Workbench
    version *below* the floor -- which would mean we declare a requirement the
    target environment does not meet.
    """
    requirements = config.get("project", {}).get("dependencies", [])
    print("\nRuntime floors (the Workbench must satisfy these)")
    print(_row("package", "declared", "workbench", ""))
    violations = 0
    for raw in requirements:
        requirement = Requirement(raw)
        name = normalize(requirement.name)
        declared = str(requirement.specifier) or "(any)"
        found = freeze.get(name)
        if found is None:
            note = "not in freeze" if name not in freeze else "unknown version"
            print(_row(name, declared, "-", note))
            continue
        if requirement.specifier.contains(found, prereleases=True):
            print(_row(name, declared, found, "ok"))
        else:
            print(_row(name, declared, found, "BELOW FLOOR"))
            violations += 1
    return violations


def check_python(config: dict, python_version: str | None) -> int:
    """Report ``requires-python`` against the Workbench interpreter."""
    declared = config.get("project", {}).get("requires-python")
    print("\nPython")
    if not declared:
        print("  requires-python is not declared")
        return 0
    if not python_version:
        print(_row("requires-python", declared, "-", "pass --python-version"))
        return 0

    requirement = Requirement(f"python{declared}")
    try:
        satisfied = requirement.specifier.contains(python_version)
    except InvalidVersion:
        print(_row("requires-python", declared, python_version, "unparseable"))
        return 1
    status = "ok" if satisfied else "WORKBENCH IS BELOW FLOOR"
    print(_row("requires-python", declared, python_version, status))
    return 0 if satisfied else 1


def report_watchlist(freeze: dict[str, str | None]) -> None:
    """Print Workbench versions for packages we deliberately leave unpinned."""
    print("\nWatchlist (not pinned; informational only)")
    for name in WATCHLIST:
        found = freeze.get(normalize(name))
        if found:
            print(f"  {name:<24} {found}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Check the Workbench pins in pyproject.toml against a pip freeze."
        )
    )
    parser.add_argument(
        "freeze", type=Path, help="path to a pip freeze from the Workbench"
    )
    parser.add_argument(
        "--python-version",
        help="Workbench Python version, e.g. 3.11.8 (pip freeze omits it)",
    )
    parser.add_argument(
        "--pyproject",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "pyproject.toml",
        help="path to pyproject.toml (defaults to the repository's own)",
    )
    args = parser.parse_args()

    if not args.freeze.is_file():
        print(f"No such freeze file: {args.freeze}", file=sys.stderr)
        return 2

    config = tomllib.loads(args.pyproject.read_text(encoding="utf-8"))
    freeze = parse_freeze(args.freeze.read_text(encoding="utf-8"))

    print(f"Workbench freeze: {args.freeze} ({len(freeze)} packages)")
    print(f"Compared against: {args.pyproject}")

    problems = check_mirror_pins(collect_mirror_pins(config), freeze)
    problems += check_runtime_floors(config, freeze)
    problems += check_python(config, args.python_version)
    report_watchlist(freeze)

    print()
    if problems:
        print(
            f"{problems} pin(s) no longer match the Workbench. Update "
            "pyproject.toml, or record why the difference is expected."
        )
        return 1
    print("Every pin still matches the Workbench.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
