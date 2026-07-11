"""
Unit tests for the general `_utils.py` submodule.

These tests use mocking to isolate the functions from dependencies on the
local file system, environment variables, and GCS.
"""
from unittest.mock import MagicMock

import pytest

from aoutools.prs._utils import (
    _stage_local_file_to_gcs,
    _standardize_chromosome_column,
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

    def test_overwrites_if_file_exists_on_gcs(self, mocker):
        """
        Tests that the copy still runs when the target already exists on GCS.
        Staging always overwrites so a stale file is never silently reused.
        """
        mocker.patch('aoutools.prs._utils.hfs.exists', return_value=True)
        local_path = '/tmp/local_file.txt'
        result = _stage_local_file_to_gcs(local_path, 'my-data')
        expected_gcs_path = 'gs://fake-bucket/data/my-data/local_file.txt'
        assert result == expected_gcs_path
        self.mock_hfs_copy.assert_called_once_with(
            f'file://{local_path}', expected_gcs_path
        )

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

    def test_annotates_column_with_prefixing_expression(self, mocker):
        """
        Tests that the column is annotated with a per-row expression that
        conditionally prepends the 'chr' prefix.

        The transform is vectorized (applied to every row via `hl.if_else`)
        rather than inferred from a single sampled value, so no `.take()`
        driver round-trip should occur.
        """
        mock_hl = mocker.patch('aoutools.prs._utils.hl', MagicMock())
        mock_table = MagicMock()

        _standardize_chromosome_column(mock_table)

        # The column is standardized in a single vectorized annotate...
        mock_table.annotate.assert_called_once()
        assert 'chr' in mock_table.annotate.call_args.kwargs
        # ...built from a per-row conditional, and the raw contig is cast to a
        # string so integer contigs are tolerated.
        mock_hl.if_else.assert_called_once()
        mock_hl.str.assert_called_once_with(mock_table.chr)
        # No sampling pass: the old `.count()`/`.take()` gate is gone.
        mock_table.count.assert_not_called()
        mock_table.select.assert_not_called()
