"""
This module defines a configuration class for Polygenic Risk Score (PRS)
calculation.
"""

from dataclasses import dataclass
from typing import Optional, Union, Sequence
import hail as hl

@dataclass
class PRSConfig:
    # pylint: disable=too-many-instance-attributes
    """
    A configuration class for Polygenic Risk Score (PRS) calculation.

    Attributes
    ----------
    chunk_size : int, default=20000
        The number of variants to include in each processing chunk.
    samples_to_keep : Union[
        hl.Table,
        Sequence[str],
        Sequence[int],
        str,
        int,
    ], optional
        A collection of sample IDs to keep. Accepts a Hail Table, or a Python
        list, set, tuple of strings or integers, or a single string or integer.
        If None, all samples are retained.
    weight_col_name : str, default='weight'
        The column name in weights table that contains effect sizes or weights.
    log_transform_weight : bool, default=False
        If True, applies a natural log transformation to the weight column.
        Useful when weights are odds ratios (OR), since PRS assumes additive
        effects on the log-odds scale.
    include_n_matched : bool, default=False
        If True, adds a column 'n_matched' with the number of variants matched
        between weights table and VDS. This option has a performance cost and
        should be used only when necessary.
    sample_id_col : str, default='person_id'
        The column name to use for sample IDs in the final output table.
    split_multi : bool, default=True
        If True, splits multi-allelic variants in VDS into bi-allelic variants
        prior to calculation.
    ref_is_effect_allele : bool, default=False
        If True, assumes effect allele in weights file corresponds to reference
        allele in VDS. Used only when `split_multi` is True.
    strict_allele_match : bool, default=True
        Used only when `split_multi` is False. If True, enforces that one
        allele in weights table matches reference allele in VDS and other
        allele is a valid alternate. If False, only effect allele is checked to
        correspond to either reference or alternate allele, and other allele is
        not verified.
    detailed_timings : bool, default=False
        If True, logs timing information for each major step. Helpful for
        diagnosing performance issues.
    """
    chunk_size: int = 20000
    samples_to_keep: Optional[
        Union[
            hl.Table,
            Sequence[str],
            Sequence[int],
            str,
            int
        ]
    ] = None
    weight_col_name: str = 'weight'
    log_transform_weight: bool = False
    include_n_matched: bool = False
    sample_id_col: str = 'person_id'
    split_multi: bool = True
    ref_is_effect_allele: bool = False
    strict_allele_match: bool = True
    detailed_timings: bool = False
