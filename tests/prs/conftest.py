"""
This module contains shared fixtures for the PRS test suite.

Fixtures defined here are automatically available to any test file within the
same directory or subdirectories, without needing to be imported.
"""

import pytest
from aoutools.prs import PRSConfig


@pytest.fixture
def default_prs_config():
    """Provides a default, reusable PRSConfig object for tests."""
    return PRSConfig()
