"""
This module defines a configuration class for Polygenic Risk Score (PRS)
calculation.
"""

# Annotations are strings, never evaluated at import. Required because Sphinx
# mocks `hail` when building the docs, and a mocked `hl.Table` does not support
# the PEP 604 `|` operator -- evaluating these annotations eagerly would break
# autodoc (real hail is unaffected).
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import hail as hl


@dataclass
class PRSConfig:
    # pylint: disable=too-many-instance-attributes
    """
    A configuration class for Polygenic Risk Score (PRS) calculation.

    Attributes
    ----------
    chunk_size : int, default 20000
        The number of variants to include in each processing chunk.
    samples_to_keep : hl.Table | Sequence[str] | Sequence[int] | str | int, \
optional
        A collection of sample IDs to keep. Accepts a Hail Table, or a Python
        list, set, tuple of strings or integers, or a single string or integer.
        If None, all samples are retained.
    weight_col_name : str, default 'weight'
        The column name in weights table that contains effect sizes or weights.
    log_transform_weight : bool, default False
        If True, applies a natural log transformation to the weight column.
        Useful when weights are odds ratios (OR), since PRS assumes additive
        effects on the log-odds scale.
    include_n_matched : bool, default False
        If True, adds a column 'n_matched' with the number of variants matched
        between weights table and VDS. This option has a performance cost and
        should be used only when necessary.
    sample_id_col : str, default 'person_id'
        The column name to use for sample IDs in the final output table.
    detailed_timings : bool, default False
        If True, adds a per-stage timing breakdown to the INFO log, useful for
        diagnosing performance issues. This is independent of the log level:
        use it to profile, and separately lower the ``aoutools`` logger to
        ``DEBUG`` if you want step-by-step detail.
    """

    chunk_size: int = 20000
    samples_to_keep: (
        hl.Table | Sequence[str] | Sequence[int] | str | int | None
    ) = None
    weight_col_name: str = "weight"
    log_transform_weight: bool = False
    include_n_matched: bool = False
    sample_id_col: str = "person_id"
    detailed_timings: bool = False
