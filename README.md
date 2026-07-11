> **Disclaimer:** This project is **not affiliated with, endorsed by, or
> sponsored by** the *All of Us Research Program*. The software is provided **as
> is**, without warranty. It is in an **early stage of development**, and its
> functions, APIs, and signatures may change periodically.

# aoutools: Tools for All of Us Researcher Workbench

aoutools is a Python library designed to simplify common analysis tasks on the
All of Us Researcher Workbench.

The initial release focuses on the `aoutools.prs` submodule, which offers
convenient functions for:

1.  **Reading PRS Weights Files:** A flexible reader that can handle various
    file formats, both with and without headers.
2.  **Calculating PRS:** A cost-efficient strategy for calculating PRS directly
    on the All of Us VDS, with support for batch mode to calculate multiple
    scores at once.

You can install `aoutools` via `pip` using either the Python Package Index
(PyPI) or its GitHub repository. On the All of Us Researcher Workbench, you can
run the following commands directly in a Jupyter Notebook cell.

```bash
# 1. From PyPI
!pip install aoutools

# 2. From Github
!pip install git+https://github.com/dokyoonkimlab/aoutools.git
```

Please check the online [aoutools
Documentation](https://aoutools.readthedocs.io) for how-to guides and API
reference.

## Development

Assumes [pixi](https://pixi.sh) and [direnv](https://direnv.net) are already
installed.

```bash
pixi install          # build the env from pixi.lock
direnv allow          # auto-activate it on cd into the repo (once per clone)
pixi run setup-hooks  # install the ruff pre-commit hooks (once per clone)
```

The test suite mocks `hail` but still imports it, so it runs only in the
linux-64 `ci` environment:

```bash
pixi run -e ci test   # full test suite
pixi run docs         # build the Sphinx HTML docs
pixi run lint         # ruff lint + format check (what CI runs)
pixi run format       # auto-fix and format
```

Formatting is enforced by `ruff format`. The bulk-reformat commit is listed in
`.git-blame-ignore-revs`; to skip it in local blame, run
`git config blame.ignoreRevsFile .git-blame-ignore-revs` (GitHub applies it
automatically).
