"""
Unit tests for the `_calculator_batch.py` submodule.

These tests use mocking to isolate the batch calculation workflow from any
real Hail/Spark or GCS dependencies.
"""
import pytest
from unittest.mock import MagicMock, patch

from aoutools.prs._calculator_batch import (
    calculate_prs_batch,
    _prepare_batch_weights_data,
    _calculate_prs_chunk_batch
)
from aoutools.prs import PRSConfig


@pytest.fixture
def mock_batch_dependencies(mocker):
    """
    A fixture to mock all external dependencies for the batch calculator tests.
    """
    # Mock entire libraries
    mocker.patch('aoutools.prs._calculator_batch.hl', MagicMock())
    mocker.patch('aoutools.prs._calculator_batch.hfs', MagicMock())
    mocker.patch('aoutools.prs._calculator_batch.pd', MagicMock())

    # Mock all imported helper functions to test orchestration
    mocker.patch('aoutools.prs._calculator_batch._prepare_samples_to_keep')
    mocker.patch(
        'aoutools.prs._calculator_batch._prepare_weights_for_chunking',
        return_value=(MagicMock(), 2)
    )
    mocker.patch(
        'aoutools.prs._calculator_batch._process_chunks_batch',
        return_value=[MagicMock()]
    )
    mocker.patch('aoutools.prs._calculator_batch._aggregate_and_export_batch')


class TestCalculatePRSBatch:
    """
    Tests for the main `calculate_prs_batch` entrypoint function.
    """

    def test_main_workflow_calls_correct_functions(
        self, mocker, mock_batch_dependencies
    ):
        """
        Tests that the main batch function calls the key helper functions
        in the correct sequence.
        """
        # Arrange: Mock the main helper functions to assert their calls
        mock_prep_batch_weights = mocker.patch(
            'aoutools.prs._calculator_batch._prepare_batch_weights_data',
            return_value=(MagicMock(), MagicMock())
        )
        mock_prep_chunking = mocker.patch(
            'aoutools.prs._calculator_batch._prepare_weights_for_chunking',
            return_value=(MagicMock(), 2)
        )
        mock_process_chunks = mocker.patch(
            'aoutools.prs._calculator_batch._process_chunks_batch',
            return_value=[MagicMock()]
        )
        mock_agg_export = mocker.patch(
            'aoutools.prs._calculator_batch._aggregate_and_export_batch'
        )

        # Act: Call the main function
        calculate_prs_batch(
            weights_tables_map={'score1': MagicMock()},
            vds=MagicMock(),
            output_path='gs://fake/path.tsv',
            config=PRSConfig()
        )

        # Assert: Check that each major step was called once
        mock_prep_batch_weights.assert_called_once()
        mock_prep_chunking.assert_called_once()
        mock_process_chunks.assert_called_once()
        mock_agg_export.assert_called_once()

    def test_aborts_if_no_variants_found(self, mocker, mock_batch_dependencies):
        """
        Tests that the function returns None and aborts if the initial
        preparation step finds no variants.
        """
        # Arrange: Make the prep function return None for loci_to_keep
        mocker.patch(
            'aoutools.prs._calculator_batch._prepare_batch_weights_data',
            return_value=(MagicMock(), None)
        )
        mock_process_chunks = mocker.patch(
            'aoutools.prs._calculator_batch._process_chunks_batch'
        )

        # Act
        result = calculate_prs_batch(
            weights_tables_map={'score1': MagicMock()},
            vds=MagicMock(),
            output_path='gs://fake/path.tsv'
        )

        # Assert
        assert result is None
        mock_process_chunks.assert_not_called()


class TestBatchHelpers:
    """
    Tests for internal helper functions in the batch calculation workflow.
    """

    def test_prepare_batch_weights_data(self, mocker):
        """
        Tests that the batch preparation function processes each table and
        creates a union of loci.
        """
        # Arrange
        mock_hl = mocker.patch('aoutools.prs._calculator_batch.hl', MagicMock())
        mock_validate = mocker.patch(
            'aoutools.prs._calculator_batch._validate_and_prepare_weights_table',
            side_effect=lambda weights_table, config: weights_table
        )
        mock_orient = mocker.patch(
            'aoutools.prs._calculator_batch._orient_weights_for_split',
            side_effect=lambda ht, config: ht
        )

        # Act
        weights_map = {'s1': MagicMock(), 's2': MagicMock()}
        _prepare_batch_weights_data(weights_map, PRSConfig(split_multi=True))

        # Assert
        assert mock_validate.call_count == 2
        assert mock_orient.call_count == 2
        mock_hl.Table.union.assert_called_once()

    def test_calculate_prs_chunk_batch(self, mocker):
        """
        Tests that the batch chunk calculator builds annotations and
        aggregators for each score.
        """
        # Arrange
        mocker.patch('aoutools.prs._calculator_batch.hl', MagicMock())
        mock_build_rows = mocker.patch(
            'aoutools.prs._calculator_batch._build_row_annotations',
            return_value={'anno1': MagicMock()}
        )
        mock_build_agg = mocker.patch(
            'aoutools.prs._calculator_batch._build_prs_agg_expr',
            return_value=MagicMock()
        )
        mock_vds = MagicMock()
        mock_mt = mock_vds.variant_data

        # Ensure annotate_rows returns the same mock object so that the
        # subsequent chained calls are made on the configured mock.
        mock_mt.annotate_rows.return_value = mock_mt

        mock_mt.select_cols().cols().select_globals.return_value = "FinalTable"

        # Act
        weights_map = {'score1': MagicMock(), 'score2': MagicMock()}
        result = _calculate_prs_chunk_batch(
            vds=mock_vds,
            weights_tables_map=weights_map,
            prepared_weights=MagicMock(),
            config=PRSConfig(split_multi=False)
        )

        # Assert
        assert result == "FinalTable"
        mock_build_rows.assert_called_once()
        assert mock_build_agg.call_count == 2
