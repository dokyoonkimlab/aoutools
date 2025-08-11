"""
Unit tests for the `_config.py` submodule.

These tests verify that the PRSConfig dataclass is instantiated correctly
and that its default attributes are set as expected.
"""
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
    assert config.weight_col_name == 'weight'
    assert config.log_transform_weight is False
    assert config.include_n_matched is False
    assert config.sample_id_col == 'person_id'
    assert config.split_multi is True
    assert config.ref_is_effect_allele is False
    assert config.strict_allele_match is True
    assert config.detailed_timings is False


def test_prs_config_custom_values():
    """
    Tests that custom values can be successfully assigned during
    instantiation.
    """
    # Arrange: Create an instance with custom, non-default values
    config = PRSConfig(
        chunk_size=100,
        weight_col_name='BETA',
        split_multi=False,
        detailed_timings=True
    )

    # Assert: Check that the custom values were assigned correctly
    assert config.chunk_size == 100
    assert config.weight_col_name == 'BETA'
    assert config.split_multi is False
    assert config.detailed_timings is True
