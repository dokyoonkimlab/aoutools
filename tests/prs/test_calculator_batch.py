"""
Unit tests for the `_calculator_batch.py` submodule.

These tests use mocking to isolate the batch calculation workflow from any
real Hail/Spark or GCS dependencies.
"""

from unittest.mock import MagicMock

import pytest

from aoutools.prs import PRSConfig
from aoutools.prs._calculator_batch import (
    _calculate_prs_chunk_batch,
    _prepare_batch_weights_data,
    calculate_prs_batch,
)


@pytest.fixture
def mock_batch_dependencies(mocker):
    """
    A fixture to mock all external dependencies for the batch calculator tests.
    """
    # Mock entire libraries
    mocker.patch("aoutools.prs._calculator_batch.hl", MagicMock())
    mocker.patch("aoutools.prs._calculator_batch.hfs", MagicMock())
    mocker.patch("aoutools.prs._calculator_batch.pd", MagicMock())

    # Mock all imported helper functions to test orchestration
    mocker.patch("aoutools.prs._calculator_batch._prepare_samples_to_keep")
    mocker.patch(
        "aoutools.prs._calculator_batch._prepare_weights_for_chunking",
        return_value=(MagicMock(), 2),
    )
    mocker.patch(
        "aoutools.prs._calculator_batch._process_chunks_batch",
        return_value=[MagicMock()],
    )
    mocker.patch("aoutools.prs._calculator_batch._aggregate_and_export_batch")


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
            "aoutools.prs._calculator_batch._prepare_batch_weights_data",
            return_value=(MagicMock(), MagicMock()),
        )
        mock_prep_chunking = mocker.patch(
            "aoutools.prs._calculator_batch._prepare_weights_for_chunking",
            return_value=(MagicMock(), 2),
        )
        mock_process_chunks = mocker.patch(
            "aoutools.prs._calculator_batch._process_chunks_batch",
            return_value=[MagicMock()],
        )
        mock_agg_export = mocker.patch(
            "aoutools.prs._calculator_batch._aggregate_and_export_batch"
        )

        # Act: Call the main function
        calculate_prs_batch(
            weights_tables_map={"score1": MagicMock()},
            vds=MagicMock(),
            output_path="gs://fake/path.tsv",
            config=PRSConfig(),
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
            "aoutools.prs._calculator_batch._prepare_batch_weights_data",
            return_value=(MagicMock(), None),
        )
        mock_process_chunks = mocker.patch(
            "aoutools.prs._calculator_batch._process_chunks_batch"
        )

        # Act
        result = calculate_prs_batch(
            weights_tables_map={"score1": MagicMock()},
            vds=MagicMock(),
            output_path="gs://fake/path.tsv",
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
        mock_hl = mocker.patch("aoutools.prs._calculator_batch.hl", MagicMock())
        mock_validate = mocker.patch(
            "aoutools.prs._calculator_batch._validate_and_prepare_weights_table",
            side_effect=lambda weights_table, config: weights_table,
        )
        mock_group = mocker.patch(
            "aoutools.prs._calculator_batch._group_weights_by_locus",
            side_effect=lambda ht: ht,
        )

        # Act
        weights_map = {"s1": MagicMock(), "s2": MagicMock()}
        _prepare_batch_weights_data(weights_map, PRSConfig())

        # Assert
        assert mock_validate.call_count == 2
        assert mock_group.call_count == 2
        mock_hl.Table.union.assert_called_once()

    def test_calculate_prs_chunk_batch(self, mocker):
        """
        Tests that the batch chunk calculator builds annotations, per-score
        row offsets, and per-score entry aggregators.
        """
        # Arrange
        mocker.patch("aoutools.prs._calculator_batch.hl", MagicMock())
        mock_split = mocker.patch(
            "aoutools.prs._calculator_batch._split_multi_with_total_dosage",
            return_value=MagicMock(),
        )
        mock_build_rows = mocker.patch(
            "aoutools.prs._calculator_batch._build_row_annotations",
            return_value={"anno1": MagicMock()},
        )
        mock_build_entry = mocker.patch(
            "aoutools.prs._calculator_batch._build_prs_entry_term",
            return_value=MagicMock(),
        )
        mock_build_offset = mocker.patch(
            "aoutools.prs._calculator_batch._build_row_offset_expr",
            return_value=MagicMock(),
        )
        mock_vds = MagicMock()
        # The MatrixTable always comes from the split VDS, which also carries
        # each entry's pre-split total non-ref count through the split.
        mock_mt = mock_split.return_value

        # Ensure annotate_rows and persist return the same mock object so that
        # the subsequent chained calls are made on the configured mock. The
        # chunk is persisted because the offset pass runs alongside the scores.
        mock_mt.annotate_rows.return_value = mock_mt
        mock_mt.persist.return_value = mock_mt

        # The chunk now renames the sample column and materializes to a
        # DataFrame before returning.
        mock_mt.select_cols().cols().select_globals().rename().to_pandas.return_value = (  # noqa: E501
            "FinalDF"
        )

        # Act
        weights_map = {"score1": MagicMock(), "score2": MagicMock()}
        result = _calculate_prs_chunk_batch(
            vds=mock_vds,
            weights_tables_map=weights_map,
            prepared_weights=MagicMock(),
            config=PRSConfig(),
        )

        # Assert
        assert result == "FinalDF"
        # The persisted chunk is released once its result is materialized.
        mock_mt.unpersist.assert_called_once()
        mock_build_rows.assert_called_once()
        # One entry term and one row-offset expression per score.
        assert mock_build_entry.call_count == 2
        assert mock_build_offset.call_count == 2
        # The row quantities are one rows-only aggregation, not a per-score
        # `aggregate_rows` folded into select_cols.
        mock_mt.rows().aggregate.assert_called_once()
        mock_mt.aggregate_rows.assert_not_called()

    def test_n_matched_does_not_add_an_entry_pass(self, mocker):
        """`include_n_matched=True` must not add a pass over the *entries* --
        folding a per-score `mt.aggregate_rows` into `select_cols` re-ran the
        whole chunk. The offsets and counts are instead one rows-only
        `mt.rows().aggregate`, leaving a single entry pass. See the single-score
        twin in `test_calculator.py`.
        """
        mocker.patch("aoutools.prs._calculator_batch.hl", MagicMock())
        mock_split = mocker.patch(
            "aoutools.prs._calculator_batch._split_multi_with_total_dosage",
            return_value=MagicMock(),
        )
        mocker.patch(
            "aoutools.prs._calculator_batch._build_row_annotations",
            return_value={"anno1": MagicMock()},
        )
        mocker.patch(
            "aoutools.prs._calculator_batch._build_prs_entry_term",
            return_value=MagicMock(),
        )
        mocker.patch(
            "aoutools.prs._calculator_batch._build_row_offset_expr",
            return_value=MagicMock(),
        )
        mock_mt = mock_split.return_value
        mock_mt.annotate_rows.return_value = mock_mt
        mock_mt.persist.return_value = mock_mt

        _calculate_prs_chunk_batch(
            vds=MagicMock(),
            weights_tables_map={"score1": MagicMock(), "score2": MagicMock()},
            prepared_weights=MagicMock(),
            config=PRSConfig(include_n_matched=True),
        )

        # The offsets and counts are one rows-only aggregation -- not folded
        # per-score into select_cols -- and there is no separate count pass.
        mock_mt.rows().aggregate.assert_called_once()
        mock_mt.aggregate_rows.assert_not_called()
        mock_mt.count_rows.assert_not_called()
        # The persisted chunk is released once its result is materialized.
        mock_mt.unpersist.assert_called_once()

    def test_effect_allele_is_alt_skips_the_offset_pass(self, mocker):
        """`effect_allele_is_alt=True` asserts no variant is REF-effect, so the
        per-score offsets are zero and their rows-only reduction is skipped.
        With `include_n_matched=False` (as here) no row reduction runs, so the
        per-sample scores are the only pass; `n_matched` would keep a second
        pass. The per-score entry terms are still built.
        """
        mocker.patch("aoutools.prs._calculator_batch.hl", MagicMock())
        mock_split = mocker.patch(
            "aoutools.prs._calculator_batch._split_multi_with_total_dosage",
            return_value=MagicMock(),
        )
        mocker.patch(
            "aoutools.prs._calculator_batch._build_row_annotations",
            return_value={"anno1": MagicMock()},
        )
        mock_build_entry = mocker.patch(
            "aoutools.prs._calculator_batch._build_prs_entry_term",
            return_value=MagicMock(),
        )
        mock_build_offset = mocker.patch(
            "aoutools.prs._calculator_batch._build_row_offset_expr",
            return_value=MagicMock(),
        )
        mock_mt = mock_split.return_value
        mock_mt.annotate_rows.return_value = mock_mt

        _calculate_prs_chunk_batch(
            vds=MagicMock(),
            weights_tables_map={"score1": MagicMock(), "score2": MagicMock()},
            prepared_weights=MagicMock(),
            config=PRSConfig(effect_allele_is_alt=True),
        )

        # No offsets built, no row aggregation, and no materialization: a
        # single per-sample pass.
        mock_build_offset.assert_not_called()
        mock_mt.rows().aggregate.assert_not_called()
        # Nothing persisted, so nothing is unpersisted.
        mock_mt.persist.assert_not_called()
        mock_mt.unpersist.assert_not_called()
        # The per-score entry terms are still built and aggregated once.
        assert mock_build_entry.call_count == 2
        mock_mt.select_cols().cols().select_globals.assert_called()
