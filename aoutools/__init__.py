# aoutools/__init__.py

# Defines the official version number for the package.
__version__ = "0.1.0"
__author__ = "Jaehyun Joo"

import logging

# This makes the 'prs' submodule directly accessible after importing
# 'aoutools'.
from . import prs

# Specifies which objects are imported when a user runs 'from aoutools import
# *'.
__all__ = ["prs"]

# If the application using this library doesn't configure logging, this line
# adds a "do-nothing" handler to prevent a "No handlers could be found"
# message.
logging.getLogger(__name__).addHandler(logging.NullHandler())
