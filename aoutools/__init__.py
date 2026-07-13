# aoutools/__init__.py

import logging
from importlib.metadata import PackageNotFoundError, version

# Derive the version from installed package metadata so pyproject.toml is the
# single source of truth.
try:
    __version__ = version("aoutools")
except PackageNotFoundError:  # not installed (e.g. running from a source tree)
    __version__ = "0.0.0"

__author__ = "Jaehyun Joo"

# This makes the 'prs' submodule directly accessible after importing
# 'aoutools'.
from . import prs
from ._workbench import (
    DEFAULT_VDS_PATH,
    get_google_project,
    get_vds_path,
    get_workspace_bucket,
    init_hail,
)

# Specifies which objects are imported when a user runs 'from aoutools import
# *'.
__all__ = [
    "DEFAULT_VDS_PATH",
    "get_google_project",
    "get_vds_path",
    "get_workspace_bucket",
    "init_hail",
    "prs",
]

# If the application using this library doesn't configure logging, this line
# adds a "do-nothing" handler to prevent a "No handlers could be found"
# message.
logging.getLogger(__name__).addHandler(logging.NullHandler())
