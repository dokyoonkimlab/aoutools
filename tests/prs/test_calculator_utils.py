"""
Unit tests for the `_calculator_utils.py` submodule.

These tests use mocking to isolate the functions from any real Hail/Spark
dependencies, allowing for fast and portable execution.
"""

from unittest.mock import MagicMock

import pytest

from aoutools.prs import PRSConfig
from aoutools.prs._calculator_utils import (
    _orient_weights_for_split,
    _validate_and_prepare_weights_table,
)


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
        # Mock the entire hail module (hl) to control its behavior. Patch it in
        # both modules that touch the table: `_calculator_utils` (this function)
        # and `_utils` (where `_standardize_chromosome_column` lives), otherwise
        # the latter's `hl.str(...)` would hit real hail's typecheck on a mock.
        mock_hl = mocker.patch("aoutools.prs._calculator_utils.hl", MagicMock())
        mocker.patch("aoutools.prs._utils.hl", MagicMock())

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

        # We define a "side effect" function that dynamically creates a mock
        # column with the correct .dtype when the table is accessed via
        # `table['column_name']`.
        def getitem_side_effect(column_name):
            mock_col = MagicMock()
            if column_name in ["chr", "effect_allele", "noneffect_allele"]:
                mock_col.dtype = mock_hl.tstr
            elif column_name == "pos":
                mock_col.dtype = mock_hl.tint32
            elif column_name == "weight":
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
            "chr": None,
            "pos": None,
            "effect_allele": None,
            "noneffect_allele": None,
            "weight": None,
        }
        mock_weights_table = self._get_mock_table(mocker, all_cols)

        # Act: Run the function under test
        _validate_and_prepare_weights_table(
            mock_weights_table, default_prs_config
        )

        # Assert: Verify that the expected Hail methods were called
        mock_weights_table.rename.assert_called_once()
        mock_weights_table.key_by.assert_called_with("locus")

    def test_missing_required_column_raises_error(self, mocker):
        """
        Tests that a TypeError is raised if a required column is missing.
        """
        # Arrange: Create a mock table missing 'pos' but including 'weight'
        # to ensure we test the correct validation step.
        cols_missing_pos = {"chr": None, "effect_allele": None, "weight": None}
        mock_weights_table = self._get_mock_table(mocker, cols_missing_pos)

        # Act & Assert: Check that the correct error is raised
        with pytest.raises(
            TypeError, match="Weights table is missing required column: 'pos'"
        ):
            _validate_and_prepare_weights_table(mock_weights_table, PRSConfig())


class TestOrientWeightsForSplit:
    """
    Tests for `_orient_weights_for_split`, which sets allele orientation on the
    split-multi path.

    This is the highest-risk function in the library: it decides which allele
    the weight applies to. Get it wrong and every score is still a clean,
    plausible float -- just inverted. Nothing else catches that.

    `config.ref_is_effect_allele` is a plain Python bool, so `hl.if_else` here
    is a branch selection on a real value. The mock below reproduces that
    semantic, which lets these tests assert the orientation and sign actually
    produced rather than merely that `hl.if_else` was called.
    """

    def _orient(self, mocker, *, ref_is_effect_allele):
        """Runs the function under test and returns (table, annotate kwargs)."""
        mock_hl = mocker.patch("aoutools.prs._calculator_utils.hl", MagicMock())
        mock_hl.if_else.side_effect = (
            lambda cond, if_true, if_false: if_true if cond else if_false
        )

        mock_table = MagicMock()
        mock_table.annotate.return_value = mock_table
        mock_table.key_by.return_value = mock_table

        config = PRSConfig(ref_is_effect_allele=ref_is_effect_allele)
        _orient_weights_for_split(mock_table, config)

        return mock_table, mock_table.annotate.call_args.kwargs

    def test_effect_allele_is_alt_keeps_weight_sign(self, mocker):
        """
        Default case: the effect allele is the ALT allele.

        Dosage on the split path is `GT.n_alt_alleles()`, i.e. a count of ALT
        copies, so the effect allele belongs in the ALT slot and the weight is
        used as-is.
        """
        table, kwargs = self._orient(mocker, ref_is_effect_allele=False)

        assert kwargs["alleles"] == [
            table.noneffect_allele,
            table.effect_allele,
        ]
        assert kwargs["weight"] is table.weight

    def test_ref_is_effect_allele_puts_effect_in_ref_and_negates_weight(
        self, mocker
    ):
        """
        `ref_is_effect_allele=True`: the effect allele is the REF allele.

        The effect allele moves to the REF slot, so `GT.n_alt_alleles()` now
        counts the *noneffect* allele. The weight is negated to compensate:
        the contribution becomes `-w * n_alt` instead of the true
        `w * n_effect = w * (2 - n_alt)`.

        These differ by a constant `2w` per matched variant. Rows are filtered
        identically for every sample, so that constant is shared across samples
        and rankings/standardized scores are preserved -- but absolute scores
        are offset. Dropping the negation flips the sign of the sample-varying
        term, which inverts the score (highest genetic risk scores lowest)
        while still producing a perfectly plausible distribution.
        """
        table, kwargs = self._orient(mocker, ref_is_effect_allele=True)

        assert kwargs["alleles"] == [
            table.effect_allele,
            table.noneffect_allele,
        ]
        # `-table.weight`; identity is stable across accesses on a MagicMock.
        assert kwargs["weight"] is table.weight.__neg__.return_value
        assert kwargs["weight"] is not table.weight

    @pytest.mark.parametrize("ref_is_effect_allele", [True, False])
    def test_keys_by_locus_and_alleles(self, mocker, ref_is_effect_allele):
        """
        The split path joins the weights to the MatrixTable on
        (locus, alleles); keying on anything else silently drops every variant.
        """
        table, _ = self._orient(
            mocker, ref_is_effect_allele=ref_is_effect_allele
        )

        table.key_by.assert_called_once_with("locus", "alleles")
