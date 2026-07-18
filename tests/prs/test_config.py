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
    assert config.detailed_timings is False
    # Off by default: the offset is computed, so scores are exact for any file.
    assert config.effect_allele_is_alt is False


def test_prs_config_custom_values():
    """
    Tests that custom values can be successfully assigned during
    instantiation.
    """
    # Arrange: Create an instance with custom, non-default values
    config = PRSConfig(
        chunk_size=100,
        weight_col_name="BETA",
        include_n_matched=True,
        detailed_timings=True,
    )

    # Assert: Check that the custom values were assigned correctly
    assert config.chunk_size == 100
    assert config.weight_col_name == "BETA"
    assert config.include_n_matched is True
    assert config.detailed_timings is True


@pytest.mark.parametrize(
    "removed_param",
    ["split_multi", "strict_allele_match", "ref_is_effect_allele"],
)
def test_removed_allele_handling_params_are_rejected(removed_param):
    """
    Three knobs were removed, each because it silently lost data.

    `split_multi=False` selected a scoring path that dropped every
    homozygous-reference sample at a variant whose effect allele was the
    reference base -- which reordered the cohort. `strict_allele_match` tuned
    only that path. `ref_is_effect_allele` declared allele orientation for a
    whole file, but orientation is a per-row property, so every row of the
    opposite orientation was silently dropped from the join; it is now resolved
    per row against the VDS. See `TODO.md`.

    All three must be a hard error rather than a silently ignored keyword. A
    `split_multi=False` that was accepted and ignored would look to the caller
    like it had been honored.
    """
    with pytest.raises(TypeError):
        PRSConfig(**{removed_param: False})
