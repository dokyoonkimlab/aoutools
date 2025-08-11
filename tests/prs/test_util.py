"""
Unit tests for the general `_utils.py` submodule.

These tests use mocking to isolate the functions from dependencies on the
local file system, environment variables, and GCS.
"""
import pytest
from unittest.mock import MagicMock, patch

from aoutools.prs._utils import (
    _stage_local_file_to_gcs,
    _standardize_chromosome_column
)


class TestStageLocalFileToGCS:
    """Tests for the `_stage_local_file_to_gcs` function."""

    @pytest.fixture(autouse=True)
    def mock_dependencies(self, mocker):
        """A fixture to mock all external dependencies for this test class."""
        mocker.patch('aoutools.prs._utils.os.path.exists', return_value=True)
        mocker.patch('aoutools.prs._utils.os.path.abspath',
                     side_effect=lambda x: x)
        mocker.patch('aoutools.prs._utils.os.getenv',
                     return_value='gs://fake-bucket')
        mocker.patch('aoutools.prs._utils.hfs.exists', return_value=False)
        self.mock_hfs_copy = mocker.patch('aoutools.prs._utils.hfs.copy')

    def test_returns_gcs_path_directly(self):
        """
        Tests that the function immediately returns a path that already
        points to GCS without performing any operations.
        """
        gcs_path = 'gs://some-bucket/file.txt'
        result = _stage_local_file_to_gcs(gcs_path, 'subdir')
        assert result == gcs_path
        self.mock_hfs_copy.assert_not_called()

    def test_stages_local_file(self):
        """
        Tests that a local file is correctly staged to the expected GCS path.
        """
        local_path = '/tmp/local_file.txt'
        result = _stage_local_file_to_gcs(local_path, 'my-data')
        expected_gcs_path = 'gs://fake-bucket/data/my-data/local_file.txt'
        assert result == expected_gcs_path
        self.mock_hfs_copy.assert_called_once_with(
            f'file://{local_path}', expected_gcs_path
        )

    def test_skips_copy_if_file_exists_on_gcs(self, mocker):
        """
        Tests that the copy operation is skipped if the target file
        already exists on GCS.
        """
        # Override the default mock for hfs.exists for this specific test
        mocker.patch('aoutools.prs._utils.hfs.exists', return_value=True)
        _stage_local_file_to_gcs('/tmp/local_file.txt', 'my-data')
        self.mock_hfs_copy.assert_not_called()

    def test_raises_error_if_local_file_not_found(self, mocker):
        """
        Tests that a FileNotFoundError is raised for a non-existent
        local file.
        """
        # Override the default mock for os.path.exists for this test
        mocker.patch('aoutools.prs._utils.os.path.exists', return_value=False)
        with pytest.raises(FileNotFoundError):
            _stage_local_file_to_gcs('/tmp/non_existent.txt', 'my-data')

    def test_raises_error_if_workspace_bucket_not_set(self, mocker):
        """
        Tests that an EnvironmentError is raised if the WORKSPACE_BUCKET
        environment variable is not set.
        """
        # Override the default mock for os.getenv for this test
        mocker.patch('aoutools.prs._utils.os.getenv', return_value=None)
        with pytest.raises(EnvironmentError):
            _stage_local_file_to_gcs('/tmp/local_file.txt', 'my-data')


class TestStandardizeChromosomeColumn:
    """Tests for the `_standardize_chromosome_column` function."""

    def test_adds_chr_prefix(self, mocker):
        """
        Tests that the 'chr' prefix is correctly added when missing.
        """
        mock_hl = mocker.patch('aoutools.prs._utils.hl', MagicMock())
        mock_table = MagicMock()
        mock_table.count.return_value = 1
        # Simulate a table where the chromosome is '1' (a string)
        mock_table.select().take.return_value = [MagicMock(chr='1')]

        _standardize_chromosome_column(mock_table)

        # Assert that the table was annotated to add the prefix
        mock_table.annotate.assert_called_once()
        # Check that the expression to add 'chr' was constructed
        assert mock_hl.str.call_args[0][0] == 'chr'

    def test_skips_annotation_if_prefix_exists(self, mocker):
        """
        Tests that no annotation occurs if the 'chr' prefix already exists.
        """
        mocker.patch('aoutools.prs._utils.hl', MagicMock())
        mock_table = MagicMock()
        mock_table.count.return_value = 1
        # Simulate a table where the chromosome is already 'chr1'
        mock_table.select().take.return_value = [MagicMock(chr='chr1')]

        _standardize_chromosome_column(mock_table)

        # Assert that annotate was NOT called
        mock_table.annotate.assert_not_called()

    def test_handles_empty_table(self, mocker):
        """
        Tests that the function returns immediately for an empty table.
        """
        mocker.patch('aoutools.prs._utils.hl', MagicMock())
        mock_table = MagicMock()
        mock_table.count.return_value = 0  # Simulate an empty table

        result = _standardize_chromosome_column(mock_table)

        # Assert that the original table object is returned
        assert result is mock_table
        # Assert that no further processing was attempted
        mock_table.select.assert_not_called()
