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
    Tests for the internal `_calculate_prs_chunk` function.
    """

    def test_calculate_prs_chunk_prepares_and_aggregates(self, mocker):
        """
        Tests that the chunk calculator prepares the MatrixTable by splitting
        multi-allelics and then aggregates it. There is only one scoring path;
        the non-split alternative was removed (see `TODO.md`).
        """
        # Arrange: Mock the internal preparation function
        mocker.patch("aoutools.prs._calculator.hl", MagicMock())
        mock_prepare_split = mocker.patch(
            "aoutools.prs._calculator._prepare_mt_split"
        )
        mock_mt = MagicMock()
        mock_prepare_split.return_value = mock_mt
        # The chunk is persisted because the offset pass runs alongside the
        # score; keep the chained calls on the same mock.
        mock_mt.persist.return_value = mock_mt

        # Act
        _calculate_prs_chunk(
            weights_table=MagicMock(), vds=MagicMock(), config=PRSConfig()
        )

        # Assert
        mock_prepare_split.assert_called_once()
        # The chunk is materialized once so the offset and score do not each
        # re-run split_multi and the join over the unpersisted MT.
        mock_mt.persist.assert_called_once()
        # ...and released once the result is materialized, so the cache does
        # not accumulate across chunks.
        mock_mt.unpersist.assert_called_once()
        # The hom-ref offset is aggregated over *rows*, not entries -- the only
        # way to reach samples with no entry. Losing this would zero every
        # homozygous-reference sample at a REF-effect variant.
        mock_mt.rows().aggregate.assert_called_once()
        # Verify that the single materializing entry aggregation was called.
        mock_mt.select_cols().cols.assert_called_once()

    def test_n_matched_does_not_add_an_entry_pass(self, mocker):
        """`include_n_matched=True` must not add a pass over the *entries*.

        A standalone `mt.count_rows()`, or folding the count into `select_cols`
        as an entry-scoped `mt.aggregate_rows(...)`, re-runs the whole chunk.
        The offset and the matched count are instead reduced together over
        `mt.rows()`, a rows-only scan. This pins that: `count_rows` and
        `aggregate_rows` are never called, there is one `mt.rows().aggregate`,
        and one materializing `select_cols().cols()`.
        """
        mocker.patch("aoutools.prs._calculator.hl", MagicMock())
        mock_mt = MagicMock()
        mocker.patch(
            "aoutools.prs._calculator._prepare_mt_split", return_value=mock_mt
        )
        mock_mt.persist.return_value = mock_mt

        _calculate_prs_chunk(
            weights_table=MagicMock(),
            vds=MagicMock(),
            config=PRSConfig(include_n_matched=True),
        )

        # A separate entry pass for the count is exactly what is avoided.
        mock_mt.count_rows.assert_not_called()
        # The offset and count are a single rows-only aggregation, not two
        # entry-scoped `aggregate_rows` folded into select_cols.
        mock_mt.rows().aggregate.assert_called_once()
        mock_mt.aggregate_rows.assert_not_called()
        # Still one materializing action, not two.
        mock_mt.select_cols().cols.assert_called_once()
        # The persisted chunk is released once its result is materialized.
        mock_mt.unpersist.assert_called_once()

    def test_effect_allele_is_alt_skips_the_offset_pass(self, mocker):
        """`effect_allele_is_alt=True` asserts no variant is REF-effect, so the
        offset is zero and its rows-only reduction is skipped. With
        `include_n_matched=False` (as here) no row reduction runs at all, so the
        per-sample score is the only pass. `n_matched` would need its own row
        reduction, and would keep the second pass -- that is covered separately.
        """
        mocker.patch("aoutools.prs._calculator.hl", MagicMock())
        mock_mt = MagicMock()
        mocker.patch(
            "aoutools.prs._calculator._prepare_mt_split", return_value=mock_mt
        )

        _calculate_prs_chunk(
            weights_table=MagicMock(),
            vds=MagicMock(),
            config=PRSConfig(effect_allele_is_alt=True),
        )

        # No offset pass: the rows reduction is skipped.
        mock_mt.rows().aggregate.assert_not_called()
        mock_mt.aggregate_rows.assert_not_called()
        # Single pass, so no materialization is needed either -- and with
        # nothing persisted, nothing is unpersisted.
        mock_mt.persist.assert_not_called()
        mock_mt.unpersist.assert_not_called()
        # Only the single per-sample entry aggregation remains.
        mock_mt.select_cols().cols.assert_called_once()
