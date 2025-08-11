"""
Unit tests for the `_calculator_utils.py` submodule.

These tests use mocking to isolate the functions from any real Hail/Spark
dependencies, allowing for fast and portable execution.
"""
import pytest
from unittest.mock import MagicMock
from aoutools.prs._calculator_utils import (
    _validate_and_prepare_weights_table
)
from aoutools.prs import PRSConfig


class TestValidateAndPrepareWeightsTable:
    """
    Tests for the `_validate_and_prepare_weights_table` function.
    """

    def _get_mock_table(self, mocker, columns_in_row):
        """
        Helper function to create a fully configured mock Hail Table.

        This handles the complex setup required to mock chained method
        calls and dynamic property access (like `.dtype`).
        """
        # Mock the entire hail module (hl) to control its behavior
        mock_hl = mocker.patch('aoutools.prs._calculator_utils.hl', MagicMock())

        # Create the main mock object for the Hail Table
        mock_table = MagicMock()
        mock_table.row = columns_in_row

        # Ensure that chained calls return the same configured mock object
        mock_table.rename.return_value = mock_table
        mock_table.annotate.return_value = mock_table
        mock_table.key_by.return_value = mock_table
        mock_table.select.return_value = mock_table

        # Configure simple method return values
        mock_table.count.return_value = 1

        # By returning a value that already starts with 'chr', we prevent
        # the function from entering the `if not str(sample_chr).startswith('chr')`
        # block, which contains the problematic Hail expression.
        mock_table.select('chr').take.return_value = [MagicMock(chr='chr1')]

        # We define a "side effect" function that dynamically creates a mock
        # column with the correct .dtype when the table is accessed via
        # `table['column_name']`.
        def getitem_side_effect(column_name):
            mock_col = MagicMock()
            if column_name in ['chr', 'effect_allele', 'noneffect_allele']:
                mock_col.dtype = mock_hl.tstr
            elif column_name == 'pos':
                mock_col.dtype = mock_hl.tint32
            elif column_name == 'weight':
                mock_col.dtype = mock_hl.tfloat64
            return mock_col

        # Assign the side effect to the mock's __getitem__ method
        mock_table.__getitem__.side_effect = getitem_side_effect
        return mock_table

    def test_success_path_with_mocks(self, mocker, default_prs_config):
        """
        Tests the successful execution path of the function.
        """
        # Arrange: Use the helper to create a mock table with all columns
        all_cols = {
            'chr': None, 'pos': None, 'effect_allele': None,
            'noneffect_allele': None, 'weight': None
        }
        mock_weights_table = self._get_mock_table(mocker, all_cols)

        # Act: Run the function under test
        _validate_and_prepare_weights_table(
            mock_weights_table, default_prs_config
        )

        # Assert: Verify that the expected Hail methods were called
        mock_weights_table.rename.assert_called_once()
        mock_weights_table.key_by.assert_called_with('locus')

    def test_missing_required_column_raises_error(self, mocker):
        """
        Tests that a TypeError is raised if a required column is missing.
        """
        # Arrange: Create a mock table missing 'pos' but including 'weight'
        # to ensure we test the correct validation step.
        cols_missing_pos = {'chr': None, 'effect_allele': None, 'weight': None}
        mock_weights_table = self._get_mock_table(mocker, cols_missing_pos)

        # Act & Assert: Check that the correct error is raised
        with pytest.raises(
            TypeError, match="Weights table is missing required column: 'pos'"
        ):
            _validate_and_prepare_weights_table(
                mock_weights_table, PRSConfig()
            )
