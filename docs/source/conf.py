# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import importlib.metadata

project = "aoutools"
copyright = "2025, Jaehyun Joo"
author = "Jaehyun Joo"

# Read the version from installed package metadata, mirroring
# aoutools/__init__.py, so pyproject.toml stays the single source of truth. This
# only reads metadata -- it does not import aoutools, which would fail here
# because hail is mocked for autodoc rather than installed.
#
# Import the module, not the `version` function: Sphinx reads module-level names
# in conf.py as config values, and a bare `version` name would shadow its own
# `version` config setting with a function object.
try:
    release = importlib.metadata.version("aoutools")
except importlib.metadata.PackageNotFoundError:
    release = "0.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import os
import sys

sys.path.insert(0, os.path.abspath("../.."))
# Mock heavy imports for autodoc to avoid installing large/system-specific
# packages
autodoc_mock_imports = [
    "hail",
    "hailtop",
    "google.cloud",  # if you don't want to install GCP libs in RTD
    "pgscatalog.core",  # could also be mocked if you don't want to install it
]
extensions = [
    "sphinx.ext.autodoc",  # Pulls documentation from docstrings.
    "sphinx.ext.napoleon",  # Understands Google-style docstrings.
    "sphinx_autodoc_typehints",  # Renders type hints nicely.
    "sphinx.ext.viewcode",  # Adds links to your source code.
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = [
    "custom.css",
]
html_title = "aoutools"
