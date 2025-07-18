# aoutools/prs/__init__.py

# Import the main reader functions from the internal _reader.py module
from ._reader import (
    read_prs_weights,
    read_prscs,
)

# Import the main calculator functions from the internal _calculator.py module
from ._calculator import (
    calculate_prs,
)

# Import the main calculator functions from the internal _calculator.py module
# from ._calculator_split_multi import (
#     calculate_prs_split,
#     calculate_prs_split_batch,
# )

# This defines what a user gets when they type 'from aoutools.prs import *'
__all__ = [
    'read_prs_weights',
    'read_prscs',
    'calculate_prs',
]
