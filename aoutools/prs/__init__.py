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

# Import the main batch calculator functions from the internal _calculator_batch.py module
from ._calculator_batch import (
    calculate_prs_batch,
)

# Import the PRSConfig from the internal _config.py module
from ._config import (
    PRSConfig,
)

# Import the PGS catalog downloader from the internal _downloader.py module
from ._downloader import (
    download_pgs,
)

# This defines what a user gets when they type 'from aoutools.prs import *'
__all__ = [
    'read_prs_weights',
    'read_prscs',
    'calculate_prs',
    'calculate_prs_batch',
    'PRSConfig',
    'download_pgs',
]
