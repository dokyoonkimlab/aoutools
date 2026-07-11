"""
Unit tests for the `_config.py` submodule.

These tests verify that the PRSConfig dataclass is instantiated correctly
and that its default attributes are set as expected.
"""

import pytest

from aoutools.prs import PRSConfig


def test_prs_config_defaults():
    """
    Tests the default attribute values of the PRSConfig dataclass.
    """
    # Arrange: Create an instance of the class with default settings
    config = PRSConfig()

    # Assert: Check that each default value is correct
    assert config.chunk_size == 20000
    assert config.samples_to_keep is None
    assert config.weight_col_name == "weight"
    assert config.log_transform_weight is False
    assert config.include_n_matched is False
    assert config.sample_id_col == "person_id"
    assert config.ref_is_effect_allele is False
    assert config.detailed_timings is False


def test_prs_config_custom_values():
    """
    Tests that custom values can be successfully assigned during
    instantiation.
    """
    # Arrange: Create an instance with custom, non-default values
    config = PRSConfig(
        chunk_size=100,
        weight_col_name="BETA",
        ref_is_effect_allele=True,
        detailed_timings=True,
    )

    # Assert: Check that the custom values were assigned correctly
    assert config.chunk_size == 100
    assert config.weight_col_name == "BETA"
    assert config.ref_is_effect_allele is True
    assert config.detailed_timings is True


@pytest.mark.parametrize(
    "removed_param", ["split_multi", "strict_allele_match"]
)
def test_removed_non_split_params_are_rejected(removed_param):
    """
    The non-split scoring path was removed: it silently dropped every
    homozygous-reference sample at a variant whose effect allele was the
    reference base, which reordered the cohort. Confirmed on the real All of Us
    VDS; see `TODO.md` (Findings 1-3).

    Both knobs that selected or tuned it must now be a hard error, not a
    silently ignored keyword -- an ignored `split_multi=False` would look like
    it had been honored.
    """
    with pytest.raises(TypeError):
        PRSConfig(**{removed_param: False})
