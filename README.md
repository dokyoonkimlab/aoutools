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
    file formats, with and without headers.
2.  **Calculating PRS:** A cost-efficient strategy for calculating PRS directly
    on the All of Us VDS, with support for batch mode to calculate multiple
    scores at once.

You can install aoutools via pip on All of Us Workbench:

```bash
# Use !pip if running in a Jupyter notebook

# 1. From PyPI
pip install aoutools

# 2. From Github
pip install git+https://github.com/dokyoonkimlab/aoutools.git
```

Please check the online [aoutools
Documentation](https://aoutools.readthedocs.io) for how-to guides and API
reference.
