"""
Unit tests for the `_reader.py` submodule.

These tests use mocking to isolate the file reading and validation logic
from any real Hail/Spark or GCS dependencies.
"""

import logging
from unittest.mock import MagicMock

import pytest

from aoutools.prs._reader import (
    _check_duplicated_ids,
    _validate_alleles,
    read_prs_weights,
    read_prscs,
)


@pytest.fixture
def mock_dependencies(mocker):
    """A fixture to mock all external dependencies for the reader tests."""
    # Mock the entire hail library
    mock_hl = mocker.patch("aoutools.prs._reader.hl", MagicMock())

    # Mock the hailtop file system
    mocker.patch("aoutools.prs._reader.hfs", MagicMock())

    # Mock the utility functions imported from other modules
    mocker.patch(
        "aoutools.prs._reader._stage_local_file_to_gcs",
        side_effect=lambda x, sub_dir: x,
    )
    mocker.patch(
        "aoutools.prs._reader._standardize_chromosome_column",
        side_effect=lambda x: x,
    )
    mocker.patch(
        "aoutools.prs._reader._process_prs_weights_table",
        side_effect=lambda table, file_path, validate_alleles: table,
    )

    return mock_hl


class TestReadPrsWeights:
    """Tests for the main `read_prs_weights` entrypoint function."""

    def test_calls_header_parser(self, mocker, mock_dependencies):
        """
        Tests that the header-based parser is called when header=True.
        """
        mock_header_parser = mocker.patch(
            "aoutools.prs._reader._read_prs_weights_header"
        )
        column_map = {
            "chr": "CHR",
            "pos": "POS",
            "effect_allele": "EA",
            "noneffect_allele": "NEA",
            "weight": "BETA",
        }
        read_prs_weights(
            file_path="gs://fake/file.txt", header=True, column_map=column_map
        )
        mock_header_parser.assert_called_once()

    def test_calls_noheader_parser(self, mocker, mock_dependencies):
        """
        Tests that the index-based parser is called when header=False.
        """
        mock_noheader_parser = mocker.patch(
            "aoutools.prs._reader._read_prs_weights_noheader"
        )
        column_map = {
            "chr": 1,
            "pos": 2,
            "effect_allele": 3,
            "noneffect_allele": 4,
            "weight": 5,
        }
        read_prs_weights(
            file_path="gs://fake/file.txt", header=False, column_map=column_map
        )
        mock_noheader_parser.assert_called_once()

    def test_raises_error_on_missing_keys_in_map(self, mock_dependencies):
        """
        Tests that a ValueError is raised if column_map is missing
        required keys.
        """
        bad_map = {"chr": 1, "pos": 2}  # Missing effect_allele, etc.
        with pytest.raises(ValueError, match="missing required keys"):
            read_prs_weights("gs://f/f.t", header=False, column_map=bad_map)

    def test_raises_error_on_mismatched_map_type(self, mock_dependencies):
        """
        Tests that a TypeError is raised if column_map values have the
        wrong type for the header setting.
        """
        # header=True expects string values, but gets an int
        bad_map = {
            "chr": "CHR",
            "pos": "POS",
            "effect_allele": "EA",
            "noneffect_allele": "NEA",
            "weight": 5,  # Should be a string
        }
        with pytest.raises(TypeError, match="must be strs"):
            read_prs_weights("gs://f/f.t", header=True, column_map=bad_map)


class TestReadPrscs:
    """The deprecated PRS-CS wrapper still works, but warns and defers to
    `read_prs_weights` with the fixed header-less column layout."""

    def test_warns_and_delegates_with_the_prscs_layout(
        self, mocker, mock_dependencies
    ):
        delegate = mocker.patch("aoutools.prs._reader.read_prs_weights")

        with pytest.warns(DeprecationWarning, match="read_prs_weights"):
            read_prscs("gs://fake/prscs.txt", validate_alleles=True)

        delegate.assert_called_once_with(
            file_path="gs://fake/prscs.txt",
            header=False,
            column_map={
                "chr": 1,
                "pos": 3,
                "effect_allele": 4,
                "noneffect_allele": 5,
                "weight": 6,
            },
            delimiter="\t",
            validate_alleles=True,
        )


class TestInternalValidators:
    """Tests for internal helper functions like allele and duplicate checks."""

    def test_validate_alleles_removes_bad_rows(self, mock_dependencies, caplog):
        """
        Tests that `_validate_alleles` filters rows and WARNs about the count.

        The removed count must reach the user at WARNING level -- silently
        dropping variants with bad alleles is the failure this guards against.
        """
        mock_table = MagicMock()
        # Simulate a table with 10 rows initially, 8 after filtering
        mock_table.count.side_effect = [10, 8]
        # Make the .filter() call return the same mock table
        mock_table.filter.return_value = mock_table

        with caplog.at_level(logging.WARNING, logger="aoutools.prs._reader"):
            result_table = _validate_alleles(mock_table)

        assert result_table.filter.call_count == 1
        assert "Removed 2 variants with invalid alleles" in caplog.text
        assert any(r.levelno == logging.WARNING for r in caplog.records)

    def test_check_duplicated_ids_raises_error(self, mock_dependencies):
        """
        Tests that `_check_duplicated_ids` raises a ValueError when
        duplicates are found.
        """
        mock_table = MagicMock()
        # This is the table after .annotate().group_by().aggregate()
        mock_agg_table = mock_table.annotate().group_by().aggregate()

        # Configure the 'n' attribute on the mock aggregate table. We replace
        # it with a new mock that has its greater-than (`__gt__`) method
        # defined. This prevents the TypeError.
        mock_n_attribute = MagicMock()
        mock_n_attribute.__gt__.return_value = "mock expression for n > 1"
        mock_agg_table.n = mock_n_attribute

        # This is the table after filtering for duplicates
        mock_duplicates_table = MagicMock()
        mock_duplicates_table.count.return_value = 1  # > 0 means duplicates
        mock_duplicates_table.take.return_value = [MagicMock(variant_id="dup1")]
        mock_agg_table.filter.return_value = mock_duplicates_table

        with pytest.raises(ValueError, match="Duplicate variants found"):
            _check_duplicated_ids(mock_table)

    def test_check_duplicated_ids_passes(self, mock_dependencies, caplog):
        """
        Tests that `_check_duplicated_ids` does not raise an error when
        no duplicates are found, and narrates only at DEBUG.

        The "no duplicates" path is routine narration, so it must stay at DEBUG
        -- default INFO output is milestones only (logging policy in AGENTS.md).
        """
        mock_table = MagicMock()
        mock_agg_table = mock_table.annotate().group_by().aggregate()

        # Apply the same fix here
        mock_n_attribute = MagicMock()
        mock_n_attribute.__gt__.return_value = "mock expression for n > 1"
        mock_agg_table.n = mock_n_attribute

        # This is the table after filtering, which should be empty
        mock_no_duplicates_table = MagicMock()
        mock_no_duplicates_table.count.return_value = 0  # 0 means no duplicates
        mock_agg_table.filter.return_value = mock_no_duplicates_table

        with caplog.at_level(logging.DEBUG, logger="aoutools.prs._reader"):
            # This should run without raising an exception
            _check_duplicated_ids(mock_table)

        assert "No duplicate variants found." in caplog.text
        assert caplog.records, "expected the routine narration to be emitted"
        assert all(r.levelno == logging.DEBUG for r in caplog.records)
