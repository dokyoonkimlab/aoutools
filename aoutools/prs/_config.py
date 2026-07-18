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
    effect_allele_is_alt : bool, default False
        Set this only if you know that in your weights file the effect allele
        is the alternate (non-reference) allele for **every** variant -- the
        usual convention of harmonized GWAS summaries. It lets the scorer skip
        the extra pass over each chunk that credits homozygous-reference
        samples at variants whose effect allele is the reference base, which can
        roughly halve the run time. (Requesting ``include_n_matched`` needs its
        own pass, so it keeps the second pass even with this set.)

        Leave it False (the default) unless you are sure: with it False the
        score is exact for any file. If you set it True and a variant's effect
        allele turns out to be the reference after all, every score is short by
        the same constant, so absolute PRS values shift but rankings,
        percentiles, and z-scores are unaffected and the cohort is never
        reordered. To check whether your file qualifies, count how many of your
        variants are reference-effect (the ``validate_public_api_on_aou.ipynb``
        notebook does this); if that count is zero, setting this True changes
        the speed but not the scores.
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
    effect_allele_is_alt: bool = False
