"""
This module defines a configuration class for the Polygenic Risk Score (PRS)
calculation.
"""

from dataclasses import dataclass
import typing

@dataclass
class PRSConfig:
    """
    A data class to hold configuration parameters for the PRS calculation.

    Attributes
    ----------
    chunk_size : int, default 20000
        The number of variants to include in each processing chunk.
    samples_to_keep : hail.Table, list, set, tuple, int, or str, optional
        A collection (list, set, tuple) or Hail Table of sample IDs to keep.
    weight_col_name : str, default 'weight'
        The name of the column in `weights_table` that contains the effect
        weights.
    log_transform_weight : bool, default False
        If True, applies a natural log transformation to the weight column. Use
        this when weights are provided as odds ratios (OR), since the PRS model
        assumes additive effects on the log-odds scale.
    include_n_shared_loci : bool, default False
        If True, adds a column 'n_shared_loci' with the total number of loci
        shared between the weights table and the VDS. Note that this can impact
        performance.
    sample_id_col : str, default 'person_id'
        The desired name for the sample ID column in the final output table.
    split_multi : bool, default True
        If True (default), splits multi-allelic variants in the VDS into
        bi-allelic variants before calculation.
    ref_is_effect_allele : bool, default False
        If True, assumes the effect allele in the weights file is the reference
        allele. This is only used when `split_multi` is True to ensure
        correct allele ordering and weight direction.
    strict_allele_match : bool, default True
        This parameter is only used when `split_multi` is False. If True, it
        performs a robust check to ensure that one allele from the weights
        table is an exact match for the reference allele in the VDS, and that
        the other allele is a valid alternate allele at that site. If False,
        the check is skipped, and the join is based on the locus only, which is
        not recommended.
    detailed_timings : bool, default False
        If True, logs the duration of each major computational step. Useful for
        debugging performance bottlenecks.
    """
    chunk_size: int = 20000
    samples_to_keep: typing.Optional[typing.Any] = None
    weight_col_name: str = 'weight'
    log_transform_weight: bool = False
    include_n_shared_loci: bool = False
    sample_id_col: str = 'person_id'
    split_multi: bool = True
    ref_is_effect_allele: bool = False
    strict_allele_match: bool = True
    detailed_timings: bool = False
