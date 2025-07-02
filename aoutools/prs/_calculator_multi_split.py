import hail as hl
import typing
import logging
import os
from ._utils import (
    _log_timing,
    _prepare_samples_to_keep,
    _standardize_chromosome_column,
)

logger = logging.getLogger(__name__)


def _validate_and_prepare_weights_table_for_split(
    table: hl.Table,
    ref_is_noneffect: bool,
    weight_col_name: str,
    log_transform_weight: bool
) -> hl.Table:
    """
    Validates and prepares the weights table for the split-based PRS calculation.

    Ensures required columns and types, standardizes chromosome format, and
    keys the table by (locus, alleles) for joining with a split MatrixTable.

    Parameters
    ----------
    table : hail.Table
        The input weights table.
    ref_is_noneffect : bool
        Flag indicating if the non-effect allele is the reference.
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
    table = table.annotate(
        locus=hl.locus(table.chr, table.pos, reference_genome='GRCh38'),
        alleles=hl.if_else(
            ref_is_noneffect,
            [table.noneffect_allele, table.effect_allele],
            [table.effect_allele, table.noneffect_allele],
        ),
        weight=hl.if_else(ref_is_noneffect, table.weight, -table.weight)
    )
    table = table.key_by('locus', 'alleles')
    return table.select('weight')


def calculate_prs_split(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    samples_to_keep: typing.Optional[typing.Any] = None,
    weight_col_name: str = 'weight',
    log_transform_weight: bool = False,
    ref_is_noneffect: bool = True,
    include_n_shared_variants: bool = True,
    sample_id_col: str = 'person_id',
    output_path: typing.Optional[str] = None,
    overwrite_output: bool = True,
    detailed_timings: bool = True
) -> hl.Table:
    """
    Calculates a Polygenic Risk Score (PRS) using a split-multi strategy.
    
    This function implements a "filter-by-locus-then-split" strategy to
    efficiently process large-scale genomic data. This method provides a
    highly accurate count of variants used in the score.

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
        A collection of sample IDs to keep.
    weight_col_name : str, default 'weight'
        The name of the column in `weights_table` that contains the effect
        weights.
    log_transform_weight : bool, default False
        If True, applies a natural log transformation to the weight column.
    ref_is_noneffect : bool, default True
        Specifies the assumed allele orientation in the weights table.
    include_n_shared_variants : bool, default True
        If True, adds a column 'n_shared_variants' with the total number of
        variants (locus and alleles) shared between the weights table and
        the VDS. This value is the same for all samples.
    sample_id_col : str, default 'person_id'
        The desired name for the sample ID column in the final output table.
    output_path : str, optional
        If provided, the final PRS table will be exported as a tab-separated
        text file to this path.
    overwrite_output : bool, default True
        If True, the function will overwrite an existing file at `output_path`.
    detailed_timings : bool, default True
        If True, logs the duration of each major computational step.

    Returns
    -------
    hail.Table
        A Hail Table with the sample ID, the calculated 'prs', and by default,
        the 'n_shared_variants' count.
    """
    logger.info("Starting PRS calculation with SPLIT-MULTI strategy...")

    gcs_output_path = None
    if output_path:
        if not output_path.startswith('gs://'):
            gcs_output_path = _stage_local_file_to_gcs(
                output_path, sub_dir='prs_results'
            )
        else:
            gcs_output_path = output_path
        if not overwrite_output and hl.hadoop_exists(gcs_output_path):
            raise FileExistsError(
                f"Output path '{gcs_output_path}' already exists."
            )

    vds_to_process = vds
    if samples_to_keep is not None:
        with _log_timing("Filtering to specified samples", detailed_timings):
            samples_ht = _prepare_samples_to_keep(samples_to_keep)
            vds_to_process = hl.vds.filter_samples(vds, samples_ht)

    with _log_timing("Validating and preparing weights table", detailed_timings):
        weights_table = _validate_and_prepare_weights_table_for_split(
            table=weights_table,
            ref_is_noneffect=ref_is_noneffect,
            weight_col_name=weight_col_name,
            log_transform_weight=log_transform_weight
        )
    
    with _log_timing("Filtering VDS to variants in weights table", detailed_timings):
        loci_to_keep = weights_table.key_by('locus').select()
        filtered_variant_data = vds_to_process.variant_data.semi_join_rows(
            loci_to_keep
        )
        filtered_vds = hl.vds.VariantDataset(
            vds_to_process.reference_data, filtered_variant_data
        )

    with _log_timing("Splitting multi-allelic variants", detailed_timings):
        mt = hl.vds.to_dense_mt(hl.vds.split_multi(filtered_vds))

    with _log_timing("Annotating variants with weights", detailed_timings):
        mt = mt.annotate_rows(prs_info=weights_table[mt.row_key])
        mt = mt.filter_rows(hl.is_defined(mt.prs_info))

    n_shared_variants = mt.count_rows()
    logger.info(f"Found {n_shared_variants} variants in common to use for PRS.")

    with _log_timing("Aggregating scores per sample", detailed_timings):
        mt = mt.annotate_cols(
            prs_score=hl.agg.sum(
                hl.coalesce(mt.GT.n_alt_alleles(), 0) * mt.prs_info.weight
            )
        )

    prs_table = mt.cols()
    final_cols = {'prs': prs_table.prs_score}
    if include_n_shared_variants:
        final_cols['n_shared_variants'] = n_shared_variants
    prs_table = prs_table.select(**final_cols)
    prs_table = prs_table.rename({'s': sample_id_col})

    if gcs_output_path:
        if overwrite_output and hl.hadoop_exists(gcs_output_path):
            hl.hadoop_rm(gcs_output_path, recursive=True)
        prs_table.export(gcs_output_path, header=True, delimiter='\t')

    logger.info("PRS calculation (split strategy) complete.")
    return prs_table


def calculate_prs_split_batch(
    weights_map: typing.Dict[str, hl.Table],
    vds: hl.vds.VariantDataset,
    samples_to_keep: typing.Optional[typing.Any] = None,
    ref_is_noneffect: bool = True,
    sample_id_col: str = 'person_id',
    output_path: typing.Optional[str] = None,
    overwrite_output: bool = True,
    detailed_timings: bool = True
) -> hl.Table:
    """
    Calculates multiple Polygenic Risk Scores (PRS) in a single pass
    using the split-multi strategy.

    This function is highly efficient for calculating several PRS on the same
    cohort, as it filters and splits the VDS only once for all combined variants.

    Parameters
    ----------
    weights_map : dict[str, hail.Table]
        A dictionary where keys are the desired PRS names (e.g., 'CAD_prs')
        and values are the corresponding Hail Tables of weights. Each table
        must contain the required columns.
    vds : hail.vds.VariantDataset
        A Variant Dataset object containing the genetic data.
    samples_to_keep : hail.Table, list, set, tuple, int, or str, optional
        A collection of sample IDs to keep.
    ref_is_noneffect : bool, default True
        Specifies the assumed allele orientation in the weights table.
    sample_id_col : str, default 'person_id'
        The desired name for the sample ID column in the final output table.
    output_path : str, optional
        If provided, the final PRS table will be exported as a tab-separated
        text file to this path.
    overwrite_output : bool, default True
        If True, the function will overwrite an existing file at `output_path`.
    detailed_timings : bool, default True
        If True, logs the duration of each major computational step.

    Returns
    -------
    hail.Table
        A "wide" Hail Table with one column per calculated PRS and one
        column for the number of shared variants for each score.
    """
    logger.info(f"Starting batch PRS calculation for {len(weights_map)} scores...")

    gcs_output_path = None
    if output_path:
        if not output_path.startswith('gs://'):
            gcs_output_path = _stage_local_file_to_gcs(
                output_path, sub_dir='prs_results'
            )
        else:
            gcs_output_path = output_path
        if not overwrite_output and hl.hadoop_exists(gcs_output_path):
            raise FileExistsError(
                f"Output path '{gcs_output_path}' already exists."
            )

    vds_to_process = vds
    if samples_to_keep is not None:
        with _log_timing("Filtering to specified samples", detailed_timings):
            samples_ht = _prepare_samples_to_keep(samples_to_keep)
            vds_to_process = hl.vds.filter_samples(vds, samples_ht)

    prepared_weights = {}
    all_loci = []
    with _log_timing("Preparing all weights tables", detailed_timings):
        for score_name, weights_table in weights_map.items():
            prepared_table = _validate_and_prepare_weights_table_for_split(
                table=weights_table,
                ref_is_noneffect=ref_is_noneffect,
                weight_col_name='weight',
                log_transform_weight=False
            )
            prepared_weights[score_name] = prepared_table
            all_loci.append(prepared_table.key_by('locus').select())

    with _log_timing("Filtering VDS to all variants", detailed_timings):
        loci_to_keep = hl.Table.union(*all_loci).key_by('locus')
        logger.info(
            f"Found {loci_to_keep.count()} total unique loci across all scores."
        )
        filtered_variant_data = vds_to_process.variant_data.semi_join_rows(
            loci_to_keep
        )
        filtered_vds = hl.vds.VariantDataset(
            vds_to_process.reference_data, filtered_variant_data
        )

    with _log_timing("Splitting multi-allelic variants", detailed_timings):
        mt = hl.vds.to_dense_mt(hl.vds.split_multi(filtered_vds))

    with _log_timing("Annotating variants with all weights", detailed_timings):
        annotation_exprs = {
            f'prs_info_{score_name}': prepared_weights[score_name][mt.row_key]
            for score_name in weights_map
        }
        mt = mt.annotate_rows(**annotation_exprs)

    with _log_timing("Calculating all PRS scores simultaneously", detailed_timings):
        aggregators = {}
        for score_name in weights_map:
            prs_info = mt[f'prs_info_{score_name}']
            variant_filter = hl.is_defined(prs_info)
            aggregators[f'{score_name}_prs'] = hl.agg.filter(
                variant_filter,
                hl.agg.sum(
                    hl.coalesce(mt.GT.n_alt_alleles(), 0) * prs_info.weight
                )
            )
            aggregators[f'{score_name}_n_variants'] = hl.agg.count_where(
                variant_filter
            )
        mt = mt.annotate_cols(**aggregators)

    logger.info("Finalizing batch PRS table...")
    prs_table = mt.cols()
    final_select_exprs = {}
    for score_name in weights_map.keys():
        final_select_exprs[score_name] = prs_table[f'{score_name}_prs']
        final_select_exprs[f'{score_name}_n_variants'] = prs_table[
            f'{score_name}_n_variants'
        ]
    prs_table = prs_table.select(**final_select_exprs)
    prs_table = prs_table.rename({'s': sample_id_col})

    if gcs_output_path:
        if overwrite_output and hl.hadoop_exists(gcs_output_path):
            hl.hadoop_rm(gcs_output_path, recursive=True)
        prs_table.export(gcs_output_path, header=True, delimiter='\t')

    logger.info("Batch PRS calculation complete.")
    return prs_table
