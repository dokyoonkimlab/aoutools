"""
Unit tests for the `_calculator.py` submodule.

These tests use mocking to isolate the main calculation workflow from
any real Hail/Spark dependencies.
"""

from unittest.mock import MagicMock

import pytest

from aoutools.prs import PRSConfig
from aoutools.prs._calculator import _calculate_prs_chunk, calculate_prs


@pytest.fixture
def mock_calculator_dependencies(mocker):
    """
    A fixture to mock all external dependencies for the calculator tests.
    This includes Hail, Pandas, and all imported helper functions.
    """
    # Mock entire libraries
    mocker.patch("aoutools.prs._calculator.hl", MagicMock())
    mocker.patch("aoutools.prs._calculator.hfs", MagicMock())
    mocker.patch("aoutools.prs._calculator.pd", MagicMock())

    # Mock all imported helper functions to test orchestration
    mocker.patch("aoutools.prs._calculator._prepare_samples_to_keep")
    mocker.patch(
        "aoutools.prs._calculator._prepare_weights_for_chunking",
        return_value=(MagicMock(), 2),  # Returns (mock_table, n_chunks)
    )
    mocker.patch(
        "aoutools.prs._calculator._process_chunks",
        return_value=[MagicMock()],  # Returns a list of mock dataframes
    )
    mocker.patch("aoutools.prs._calculator._aggregate_and_export")


class TestCalculatePRS:
    """
    Tests for the main `calculate_prs` entrypoint function, focusing on
    orchestration and control flow.
    """

    def test_main_workflow_calls_correct_functions(
        self, mocker, mock_calculator_dependencies
    ):
        """
        Tests that the main function calls the key helper functions in the
        correct sequence.
        """
        # The fixture already mocks these. We just need to get a reference to
        # them for assertion, not patch them again.
        mock_prep_weights = mocker.patch(
            "aoutools.prs._calculator._prepare_weights_for_chunking",
            return_value=(MagicMock(), 2),
        )
        mock_process_chunks = mocker.patch(
            "aoutools.prs._calculator._process_chunks",
            return_value=[MagicMock()],
        )
        mock_agg_export = mocker.patch(
            "aoutools.prs._calculator._aggregate_and_export"
        )

        # Act: Call the main function
        calculate_prs(
            weights_table=MagicMock(),
            vds=MagicMock(),
            output_path="gs://fake/path.tsv",
            config=PRSConfig(),
        )

        # Assert: Check that each major step was called once
        mock_prep_weights.assert_called_once()
        mock_process_chunks.assert_called_once()
        mock_agg_export.assert_called_once()

    def test_raises_error_for_non_gcs_path(self, mock_calculator_dependencies):
        """
        Tests that a ValueError is raised if the output path is not a
        GCS path.
        """
        with pytest.raises(ValueError, match="must be a Google Cloud Storage"):
            calculate_prs(
                weights_table=MagicMock(),
                vds=MagicMock(),
                output_path="/local/path.tsv",  # Not a 'gs://' path
            )


class TestChunkProcessing:
    """
    Tests for the internal `_calculate_prs_chunk` function to verify
    it dispatches to the correct preparation function based on config.
    """

    def _get_mock_mt(self):
        """Helper to create a mock MatrixTable for testing aggregation."""
        mock_mt = MagicMock()
        # Configure the mock MT so that the multiplication inside
        # hl.agg.sum(mt.dosage * mt.weights_info.weight) returns a number. This
        # prevents the TypeError.
        mock_mt.dosage.__mul__.return_value = 1.0
        return mock_mt

    def test_calculate_prs_chunk_prepares_and_aggregates(self, mocker):
        """
        Tests that the chunk calculator prepares the MatrixTable by splitting
        multi-allelics and then aggregates it. There is only one scoring path;
        the non-split alternative was removed (see `TODO.md`).
        """
        # Arrange: Mock the internal preparation function
        mock_prepare_split = mocker.patch(
            "aoutools.prs._calculator._prepare_mt_split"
        )

        # Create a mock MatrixTable to be returned by the prep function
        mock_mt = self._get_mock_mt()
        mock_prepare_split.return_value = mock_mt

        # Act
        _calculate_prs_chunk(
            weights_table=MagicMock(), vds=MagicMock(), config=PRSConfig()
        )

        # Assert
        mock_prepare_split.assert_called_once()
        # Verify that the final aggregation step was called on the mock mt
        mock_mt.select_cols().cols.assert_called_once()
