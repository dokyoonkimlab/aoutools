import os
import typing
import logging
import hail as hl
from ._utils import (
    _log_timing,
    _standardize_chromosome_column,
    _prepare_samples_to_keep,
)

logger = logging.getLogger(__name__)


def _validate_and_prepare_weights_table(
    table: hl.Table,
    weight_col_name: str,
    log_transform_weight: bool
) -> hl.Table:
    """
    Validates and prepares a single weights table for PRS calculation.

    This function ensures the table has the required columns with the correct
    types, standardizes the chromosome format, handles the weight column,
    and keys the table by locus for joining with the VDS.

    Parameters
    ----------
    table : hail.Table
        The input weights table. Must contain 'chr', 'pos', 'effect_allele',
        'noneffect_allele', and a weight column.
    weight_col_name : str
        The name of the column containing the effect weights.
    log_transform_weight : bool
        If True, applies a natural log transformation to the weight column.

    Returns
    -------
    hail.Table
        A validated, standardized, and keyed Hail Table.
    """
    if weight_col_name not in table.row:
        raise TypeError(
            f"Specified weight column '{weight_col_name}' not found in table."
        )
    table = table.rename({weight_col_name: 'weight'})
    required_cols = {
        'chr': hl.tstr, 'pos': hl.tint32, 'effect_allele': hl.tstr,
        'noneffect_allele': hl.tstr, 'weight': hl.tfloat64,
    }
    for col, expected_type in required_cols.items():
        if col not in table.row:
            raise TypeError(f"Weights table is missing required column: '{col}'.")
        if table[col].dtype != expected_type:
            raise TypeError(f"Column '{col}' has incorrect type.")
    table = _standardize_chromosome_column(table)
    if log_transform_weight:
        table = table.annotate(weight=hl.log(table.weight))
    table = table.annotate(locus=hl.locus(table.chr, table.pos,
                                          reference_genome='GRCh38'))
    table = table.key_by('locus')
    return table.select('effect_allele', 'noneffect_allele', 'weight')


def _calculate_dosage(
    mt: hl.MatrixTable,
    score_name: str = ''
) -> hl.expr.Int32Expression:
    """
    Calculates the dosage of the effect allele for a specific score.

    This expression handles both global (GT) and local (LGT/LA) genotype
    encoding formats, which enables sparse storage of homozygous reference
    calls. It correctly computes dosage at multi-allelic sites.

    Parameters
    ----------
    mt : hail.MatrixTable
        MatrixTable annotated with a `weights_info` struct. For batch mode,
        the struct is named `weights_info_{score_name}`.
    score_name : str, optional
        The identifier for the PRS being calculated, used to access the
        correct weights annotation in batch mode.

    Returns
    -------
    hail.expr.Int32Expression
        An expression for the effect allele dosage.
    """
    weights_field_name = (
        f'weights_info_{score_name}' if score_name else 'weights_info'
    )
    effect_allele = mt[weights_field_name].effect_allele
    ref_is_effect = (effect_allele == mt.alleles[0])

    # Check for 'GT' field to handle different VDS versions by their
    # genotype encoding scheme.
    if 'GT' in mt.entry:
        # Global-indexed format: 'GT' contains indices that refer
        # directly to the global 'alleles' array.
        # Example: if GT is [0, 1] and mt.alleles is ['A', 'G', 'T'],
        # this expression reconstructs the sample's alleles as ['A', 'G'].
        alleles_expr = hl.or_missing(
            hl.is_defined(mt.GT),
            hl.array([mt.alleles[mt.GT[0]], mt.alleles[mt.GT[1]]])
        )
    else:
        # Local-indexed (sparse) format: 'LGT' indices refer to the 'LA'
        # (local-to-global) map, which then refers to 'alleles'.
        # Example: LGT=[0, 1], LA=[0, 2], mt.alleles=['A', 'C', 'G']
        # 1. LGT[0] is 0. LA[0] is 0. mt.alleles[0] is 'A'.
        # 2. LGT[1] is 1. LA[1] is 2. mt.alleles[2] is 'G'.
        # The reconstructed alleles are ['A', 'G'].
        alleles_expr = hl.or_missing(
            hl.is_defined(mt.LGT) & hl.is_defined(mt.LA),
            hl.array([
                mt.alleles[hl.or_else(mt.LA[mt.LGT[0]], 0)],
                mt.alleles[hl.or_else(mt.LA[mt.LGT[1]], 0)]
            ])
        )

    # The hl.case statement cleanly handles missing genotypes by assuming
    # they are homozygous reference.
    return hl.case() \
        .when(hl.is_missing(alleles_expr) & ref_is_effect, 2) \
        .when(hl.is_missing(alleles_expr) & ~ref_is_effect, 0) \
        .default(
            hl.or_else(alleles_expr, hl.empty_array(hl.tstr)).filter(
                lambda allele: allele == effect_allele
            ).length()
        )


def calculate_prs(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    samples_to_keep: typing.Optional[typing.Any] = None,
    weight_col_name: str = 'weight',
    log_transform_weight: bool = False,
    strict_allele_match: bool = True,
    include_n_shared_loci: bool = True,
    sample_id_col: str = 'person_id',
    output_path: typing.Optional[str] = None,
    overwrite_output: bool = True,
    detailed_timings: bool = True
) -> hl.Table:
    """Calculates a Polygenic Risk Score (PRS) without splitting multi-allelic
    variants.

    This function is optimized for performance on the All of Us VDS, where
    multi-allelic variants are not pre-split. It computes allele dosages
    directly from the multi-allelic VDS, avoiding the computationally
    expensive `split_multi` step. It is compatible with different VDS
    versions by handling both global ('GT') and local ('LGT'/'LA') genotype
    encoding schemes.

    To achieve this, the function joins the VDS and the weights table using
    only the genomic `locus` as a key. A critical consequence of this design
    is that it assumes there is only **one** variant entry per locus in the
    `weights_table`. If the weights table contains multiple entries for the
    same locus, Hail will use the first matching entry it encounters, which
    can lead to an incorrect PRS calculation. If you prefer additional
    robustness over speed, consider using the `calculate_prs_split` function.

    Parameters
    ----------
    weights_table : hail.Table
        A table of PRS weights. Must contain the following columns:
        - 'chr': str
        - 'pos': int32
        - 'effect_allele': str
        - 'noneffect_allele': str
        - A column for the effect weight (float64), specified by
          `weight_col_name`.
    vds : hail.vds.VariantDataset
        A Variant Dataset object containing the genetic data.
    samples_to_keep : hail.Table, list, set, tuple, int, or str, optional
        A collection of sample IDs to keep. If provided, the PRS will only be
        calculated for these samples. Numeric IDs will be converted to strings.
    weight_col_name : str, default 'weight'
        The name of the column in `weights_table` that contains the effect
        weights.
    log_transform_weight : bool, default False
        If True, applies a natural log transformation to the weight column.
        This should be used when weights are provided as odds ratios (OR).
    strict_allele_match : bool, default True
        If True, validates that one of the alleles from weight table is the
        reference and the other is a valid alternate allele in the VDS.
    include_n_shared_loci : bool, default True
        If True, adds a column 'n_shared_loci' with the total number of
        loci shared between the weights table and the VDS.
    sample_id_col : str, default 'person_id'
        The desired name for the sample ID column in the final output table.
    detailed_timings : bool, default True
        If True, logs the duration of each major computational step.

    Returns
    -------
    hail.Table
        A Hail Table with the sample ID, the calculated 'prs', and by default,
        the 'n_shared_loci' count.

    """
    logger.info("Starting PRS calculation...")

    if samples_to_keep is not None:
        with _log_timing("Filtering to specified samples", detailed_timings):
            samples_ht = _prepare_samples_to_keep(samples_to_keep)
            vds = hl.vds.filter_samples(vds, samples_ht)

    with _log_timing("Validating and preparing weights table", detailed_timings):
        weights_table = _validate_and_prepare_weights_table(
            weights_table, weight_col_name, log_transform_weight
        )

    with _log_timing("Filtering VDS to variants in weights table", detailed_timings):
        mt = vds.variant_data.semi_join_rows(weights_table)

    with _log_timing("Annotating variants with weights", detailed_timings):
        mt = mt.annotate_rows(weights_info=weights_table[mt.locus])

    if strict_allele_match:
        with _log_timing("Performing strict allele match", detailed_timings):
            initial_count = mt.count_rows()
            vds_alleles = hl.set(mt.alleles)
            ref_allele = mt.alleles[0]
            effect = mt.weights_info.effect_allele
            noneffect = mt.weights_info.noneffect_allele
            is_valid_pair = (
                ((effect == ref_allele) & vds_alleles.contains(noneffect)) |
                ((noneffect == ref_allele) & vds_alleles.contains(effect))
            ) & (effect != noneffect)
            mt = mt.filter_rows(is_valid_pair)
            n_removed = initial_count - mt.count_rows()
            if n_removed > 0:
                logger.warning(
                    f"Removed {n_removed} variants due to allele mismatch."
                )

    n_shared_loci = mt.count_rows()
    logger.info(f"Found {n_shared_loci} loci in common to use for PRS.")

    with _log_timing("Calculating per-variant dosage", detailed_timings):
        mt = mt.annotate_entries(dosage=_calculate_dosage(mt))

    with _log_timing("Aggregating scores per sample", detailed_timings):
        mt = mt.annotate_cols(
            prs_score=hl.agg.sum(mt.dosage * mt.weights_info.weight)
        )

    prs_table = mt.cols()
    final_cols = {'prs': prs_table.prs_score}
    if include_n_shared_loci:
        final_cols['n_shared_loci'] = n_shared_loci
    prs_table = prs_table.select(**final_cols)
    prs_table = prs_table.rename({'s': sample_id_col})

    logger.info("PRS calculation complete.")
    return prs_table


def calculate_prs_batch(
    weights_map: typing.Dict[str, hl.Table],
    vds: hl.vds.VariantDataset,
    samples_to_keep: typing.Optional[typing.Any] = None,
    strict_allele_match: bool = True,
    sample_id_col: str = 'person_id',
    output_path: typing.Optional[str] = None,
    overwrite_output: bool = True,
    detailed_timings: bool = True
) -> hl.Table:
    """
    Calculates multiple Polygenic Risk Scores (PRS) in a single pass.

    This function is highly efficient for calculating several PRS on the same
    cohort, as it filters the VDS only once for all combined variants.

    Note: it is optimized for performance on the All of Us VDS, where
    multi-allelic variants are not pre-split. It computes allele dosages
    directly from the multi-allelic VDS, avoiding the computationally expensive
    `split_multi` step. It is compatible with different VDS versions by
    handling both global ('GT') and local ('LGT'/'LA') genotype encoding
    schemes.

    To achieve this, the function joins the VDS and the weights table using
    only the genomic `locus` as a key. A critical consequence of this design is
    that it assumes there is only **one** variant entry per locus in the
    `weights_table`. If the weights table contains multiple entries for the
    same locus, Hail will use the first matching entry it encounters, which can
    lead to an incorrect PRS calculation. If you prefer additional robustness
    over speed, consider using the `calculate_prs_split_batch` function.

    Parameters
    ----------
    weights_map : dict[str, hail.Table]
        A dictionary where keys are the desired PRS names (e.g., 'CAD_prs')
        and values are the corresponding Hail Tables of weights. Each table
        must contain the required columns ('chr', 'pos', 'effect_allele',
        'noneffect_allele', and a weight column).
    vds : hail.vds.VariantDataset
        A Variant Dataset object containing the genetic data.
    samples_to_keep : hail.Table, list, set, tuple, int, or str, optional
        A collection of sample IDs to keep.
    strict_allele_match : bool, default True
        If True, validates that one of the alleles from weight table is the
        reference and the other is a valid alternate allele in the VDS.
    sample_id_col : str, default 'person_id'
        The desired name for the sample ID column in the final output table.
    detailed_timings : bool, default True
        If True, logs the duration of each major computational step.

    Returns
    -------
    hail.Table
        A "wide" Hail Table with one column per calculated PRS and one
        column for the number of shared loci for each score.

    """
    logger.info(f"Starting batch PRS calculation for {len(weights_map)} scores...")

    if samples_to_keep is not None:
        with _log_timing("Filtering to specified samples", detailed_timings):
            samples_ht = _prepare_samples_to_keep(samples_to_keep)
            vds = hl.vds.filter_samples(vds, samples_ht)

    prepared_weights = {}
    all_loci = []
    with _log_timing("Preparing all weights tables", detailed_timings):
        for score_name, weights_table in weights_map.items():
            prepared_table = _validate_and_prepare_weights_table(
                weights_table, 'weight', False
            )
            prepared_weights[score_name] = prepared_table
            all_loci.append(prepared_table.key_by('locus').select())

    with _log_timing("Filtering VDS to all variants", detailed_timings):
        loci_to_keep = hl.Table.union(*all_loci).key_by('locus')
        logger.info(
            f"Found {loci_to_keep.count()} total unique loci across all scores."
        )
        mt = vds.variant_data.semi_join_rows(loci_to_keep)

    with _log_timing("Annotating variants with all weights", detailed_timings):
        annotation_exprs = {
            f'weights_info_{score_name}': prepared_weights[score_name][mt.locus]
            for score_name in weights_map
        }
        mt = mt.annotate_rows(**annotation_exprs)

    if strict_allele_match:
        with _log_timing(
            "Performing strict allele match for all scores", detailed_timings
        ):
            vds_alleles = hl.set(mt.alleles)
            ref_allele = mt.alleles[0]
            new_annotations = {}
            for score_name in weights_map:
                weights_info = mt[f'weights_info_{score_name}']
                effect = weights_info.effect_allele
                noneffect = weights_info.noneffect_allele
                is_valid_pair = (
                    ((effect == ref_allele) & vds_alleles.contains(noneffect)) |
                    ((noneffect == ref_allele) & vds_alleles.contains(effect))
                ) & (effect != noneffect)
                new_annotations[f'weights_info_{score_name}'] = hl.if_else(
                    is_valid_pair,
                    weights_info,
                    hl.missing(weights_info.dtype)
                )
            mt = mt.annotate_rows(**new_annotations)

    with _log_timing("Calculating all PRS scores simultaneously", detailed_timings):
        aggregators = {}
        for score_name in weights_map:
            weights_info = mt[f'weights_info_{score_name}']
            variant_filter = hl.is_defined(weights_info)
            dosage_expr = _calculate_dosage(mt, score_name)
            aggregators[f'{score_name}_prs'] = hl.agg.filter(
                variant_filter,
                hl.agg.sum(dosage_expr * weights_info.weight)
            )
            aggregators[f'{score_name}_n_loci'] = hl.agg.count_where(
                variant_filter
            )
        mt = mt.annotate_cols(**aggregators)

    logger.info("Finalizing batch PRS table...")
    prs_table = mt.cols()
    final_select_exprs = {}
    for score_name in weights_map.keys():
        final_select_exprs[score_name] = prs_table[f'{score_name}_prs']
        final_select_exprs[f'{score_name}_n_loci'] = prs_table[
            f'{score_name}_n_loci'
        ]
    prs_table = prs_table.select(**final_select_exprs)
    prs_table = prs_table.rename({'s': sample_id_col})

    logger.info("Batch PRS calculation complete.")
    return prs_table
