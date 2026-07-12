"""
Unit tests for the `_calculator_utils.py` submodule.

These tests use mocking to isolate the functions from any real Hail/Spark
dependencies, allowing for fast and portable execution.
"""

from unittest.mock import MagicMock

import pytest

from aoutools.prs import PRSConfig
from aoutools.prs._calculator_utils import (
    _key_weights_by_variant,
    _orient_weight_and_offset,
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


class TestKeyWeightsByVariant:
    """
    Tests for `_key_weights_by_variant`, which builds the join key.

    The key is the *unordered* allele pair, so a weights row matches its variant
    whichever way round the effect allele is written. Keying on the ordered pair
    -- as the old `ref_is_effect_allele` code did -- silently drops every row of
    the opposite orientation. Keying on the locus alone would instead match the
    wrong variant when a file names two variants at one position.
    """

    def test_keys_by_locus_and_the_sorted_allele_pair(self, mocker):
        mock_hl = mocker.patch("aoutools.prs._calculator_utils.hl", MagicMock())
        mock_table = MagicMock()
        mock_table.annotate.return_value = mock_table
        mock_table.key_by.return_value = mock_table

        _key_weights_by_variant(mock_table)

        # The pair is canonicalized, not laid out as [ref, alt].
        mock_hl.sorted.assert_called_once_with(
            [mock_table.effect_allele, mock_table.noneffect_allele]
        )
        kwargs = mock_table.annotate.call_args.kwargs
        assert kwargs["alleles"] is mock_hl.sorted.return_value
        mock_table.key_by.assert_called_once_with("locus", "alleles")

    def test_does_not_touch_the_weight(self, mocker):
        """Orientation is resolved later, per row, against the VDS's own REF.
        The weight must arrive at that point unmodified -- negating it here, as
        the old code did, would double-negate."""
        mocker.patch("aoutools.prs._calculator_utils.hl", MagicMock())
        mock_table = MagicMock()
        mock_table.annotate.return_value = mock_table
        mock_table.key_by.return_value = mock_table

        _key_weights_by_variant(mock_table)

        assert "weight" not in mock_table.annotate.call_args.kwargs


class TestOrientWeightAndOffset:
    """
    Tests for `_orient_weight_and_offset`, the highest-risk function here: it
    decides which allele the weight applies to. Get the sign wrong and every
    score is still a clean, plausible float -- just inverted.

    `hl.if_else` is mocked to actually select a branch, so these assert the sign
    produced rather than merely that `hl.if_else` was called. The arithmetic
    itself is verified against real genotypes in `tests/integration/`.
    """

    def _orient(self, mocker, *, ref_is_effect):
        mock_hl = mocker.patch("aoutools.prs._calculator_utils.hl", MagicMock())
        mock_hl.if_else.side_effect = (
            lambda cond, if_true, if_false: if_true if cond else if_false
        )

        mt = MagicMock()
        weights_info = MagicMock()
        if ref_is_effect:
            # `effect_allele == alleles[0]` is True when they are the same
            # object: MagicMock's default __eq__ is identity.
            mt.alleles.__getitem__.return_value = weights_info.effect_allele

        return weights_info, _orient_weight_and_offset(mt, weights_info)

    def test_effect_allele_on_the_alt_keeps_the_weight_and_adds_no_offset(
        self, mocker
    ):
        """Dosage is `GT.n_alt_alleles()`, so an ALT effect allele needs the
        weight as-is. A hom-ref sample carries no copies of it and correctly
        contributes nothing -- hence no offset."""
        weights_info, (weight_per_alt_copy, hom_ref_offset) = self._orient(
            mocker, ref_is_effect=False
        )

        assert weight_per_alt_copy is weights_info.weight
        assert hom_ref_offset == 0.0

    def test_effect_allele_on_the_ref_negates_the_weight_and_offsets_by_2w(
        self, mocker
    ):
        """A REF effect allele means the true contribution is `w * (2 - n_alt)`,
        which is `2w - w * n_alt`. The per-entry term therefore carries the
        *negated* weight, and the `2w` becomes a row-level offset -- the only
        way to reach hom-ref samples, who have no entry to aggregate over.

        Dropping the negation inverts the sample-varying term (highest genetic
        risk scores lowest) while still producing a plausible distribution.
        Dropping the offset silently zeroes every hom-ref sample.
        """
        weights_info, (weight_per_alt_copy, hom_ref_offset) = self._orient(
            mocker, ref_is_effect=True
        )

        assert weight_per_alt_copy is weights_info.weight.__neg__.return_value
        assert weight_per_alt_copy is not weights_info.weight
        # `2.0 * weight` dispatches to the mock's __rmul__.
        assert hom_ref_offset is weights_info.weight.__rmul__.return_value
