"""
Unit tests for the `_calculator_utils.py` submodule.

These tests use mocking to isolate the functions from any real Hail/Spark
dependencies, allowing for fast and portable execution.
"""

from unittest.mock import MagicMock

import pytest

from aoutools.prs import PRSConfig
from aoutools.prs._calculator_utils import (
    _group_weights_by_locus,
    _match_weight_at_locus,
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


class TestGroupWeightsByLocus:
    """
    Tests for `_group_weights_by_locus`, which builds the join's small side.

    The weights are grouped by **locus alone** so the join to the split MT is a
    key-*prefix* join (no shuffle); alleles are matched locally afterwards by
    `_match_weight_at_locus`. Keying on the full `(locus, alleles)` instead
    would force the MT to be re-keyed to minimal alleles, shuffling every
    entry. The
    orientation-agnostic and locus-shift behaviours are exercised for real in
    `tests/integration/test_allele_matching.py`; these only pin the shape.
    """

    def _group(self, mocker):
        mock_hl = mocker.patch("aoutools.prs._calculator_utils.hl", MagicMock())
        mock_table = MagicMock()
        mock_table.filter.return_value = mock_table
        mock_grouped = MagicMock()
        mock_table.group_by.return_value = mock_grouped
        _group_weights_by_locus(mock_table)
        return mock_hl, mock_table, mock_grouped

    def test_groups_by_locus_alone(self, mocker):
        _, mock_table, _ = self._group(mocker)

        mock_table.group_by.assert_called_once_with(mock_table.locus)

    def test_collects_a_variants_array_per_locus(self, mocker):
        mock_hl, _, mock_grouped = self._group(mocker)

        kwargs = mock_grouped.aggregate.call_args.kwargs
        assert set(kwargs) == {"variants"}, (
            "each locus must carry an array of its variants for the local "
            "allele match, not a single pre-picked row"
        )
        mock_hl.agg.collect.assert_called_once()

    def test_drops_a_row_naming_the_same_allele_twice(self, mocker):
        """`A`/`A` can match no VDS row (REF != ALT), so it is filtered out
        rather than left to pollute a locus's variants array."""
        _, mock_table, _ = self._group(mocker)

        mock_table.filter.assert_called_once()

    def test_does_not_canonicalize_with_sorted(self, mocker):
        """Regression guard against the old sorted-pair key, which silently
        dropped every variant whose REF sorts after its ALT."""
        mock_hl, _, _ = self._group(mocker)

        mock_hl.sorted.assert_not_called()


class TestMatchWeightAtLocus:
    """
    Tests for `_match_weight_at_locus`, the local allele match after the join.
    """

    def test_matches_on_the_unordered_allele_set(self, mocker):
        """The match must succeed whichever way round the GWAS wrote the pair --
        orientation is resolved later, from the REF. So it compares allele
        *sets*, not ordered pairs."""
        mock_hl = mocker.patch("aoutools.prs._calculator_utils.hl", MagicMock())
        variants = MagicMock()
        canonical = MagicMock()

        _match_weight_at_locus(variants, canonical)

        # The target is the canonical alleles as a set, and the locus's
        # variants array is searched with `find` (not indexed by order).
        mock_hl.set.assert_any_call(canonical)
        variants.find.assert_called_once()


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

        weights_info = MagicMock()
        # `effect_allele == ref_allele` is True when they are the same object:
        # MagicMock's default __eq__ is identity.
        ref_allele = (
            weights_info.effect_allele if ref_is_effect else MagicMock()
        )

        return weights_info, _orient_weight_and_offset(ref_allele, weights_info)

    def test_effect_allele_on_the_alt_keeps_the_weight_and_adds_no_offset(
        self, mocker
    ):
        """Dosage is `GT.n_alt_alleles()`, so an ALT effect allele needs the
        weight as-is. A hom-ref sample carries no copies of it and correctly
        contributes nothing -- hence no offset."""
        (
            weights_info,
            (
                weight_per_alt_copy,
                hom_ref_offset,
                ref_is_effect,
            ),
        ) = self._orient(mocker, ref_is_effect=False)

        assert weight_per_alt_copy is weights_info.weight
        assert hom_ref_offset == 0.0
        assert ref_is_effect is False

    def test_effect_allele_on_the_ref_negates_the_weight_and_offsets_by_2w(
        self, mocker
    ):
        """A REF effect allele means the true contribution is
        `w * (2 - n_non_ref)`, which is `2w - w * n_non_ref`. The per-entry term
        therefore carries the *negated* weight, and the `2w` becomes a row-level
        offset -- the only way to reach hom-ref samples, who have no entry to
        aggregate over.

        Dropping the negation inverts the sample-varying term (highest genetic
        risk scores lowest) while still producing a plausible distribution.
        Dropping the offset silently zeroes every hom-ref sample.

        `ref_is_effect` is returned so `_entry_contribution` knows to count
        `n_non_ref` rather than the downcoded `GT.n_alt_alleles()`; at a
        multi-allelic site those differ. See
        `tests/integration/test_allele_matching.py`.
        """
        (
            weights_info,
            (
                weight_per_alt_copy,
                hom_ref_offset,
                ref_is_effect,
            ),
        ) = self._orient(mocker, ref_is_effect=True)

        assert weight_per_alt_copy is weights_info.weight.__neg__.return_value
        assert weight_per_alt_copy is not weights_info.weight
        # `2.0 * weight` dispatches to the mock's __rmul__.
        assert hom_ref_offset is weights_info.weight.__rmul__.return_value
        assert ref_is_effect is True
