""""PRS calculator"""

import os
import typing
import logging
from math import ceil
from uuid import uuid4
import hail as hl
import hailtop.fs as hfs
import pandas as pd
from aoutools._utils.helpers import SimpleTimer
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

    Raises
    ------
    TypeError
        If the specified `weight_col_name` does not exist, if other required
        columns are missing, or if any required column has an incorrect
        data type.
    """
    if weight_col_name not in table.row:
        raise TypeError(
            f"Specified weight column '{weight_col_name}' not found in table."
        )
    table = table.rename({weight_col_name: 'weight'})

    required_cols = {
        'chr': hl.tstr,
        'pos': hl.tint32,
        'effect_allele': hl.tstr,
        'noneffect_allele': hl.tstr,
        'weight': hl.tfloat64,
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
        locus=hl.locus(table.chr, table.pos, reference_genome='GRCh38')
    )
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
    ref_is_effect = effect_allele == mt.alleles[0]

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
        # Local-indexed format: 'LGT' indices refer to the 'LA'
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


def _prepare_mt_split(
    vds: hl.vds.VariantDataset,
    weights_table: hl.Table,
    ref_is_effect_allele: bool,
    detailed_timings: bool,
) -> hl.MatrixTable:
    """
    Prepares a MatrixTable for the split-multi PRS calculation path.

    This function takes an interval-filtered VDS, splits multi-allelic sites,
    and joins it with the weights table using a precise `(locus, alleles)`
    key. It handles allele orientation and weight direction based on the
    `ref_is_effect_allele` flag and calculates the bi-allelic dosage.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        The interval-filtered Variant Dataset.
    weights_table : hail.Table
        The chunk of the weights table.
    ref_is_effect_allele : bool
        If True, assumes the effect allele is the reference.
    detailed_timings : bool
        If True, logs the duration of computational steps.

    Returns
    -------
    hail.MatrixTable
        A prepared MatrixTable, filtered and annotated with `weights_info`
        and `dosage` for the specified effect allele.
    """
    with _log_timing(
        "Splitting multi-allelic variants and joining", detailed_timings
    ):
        split_vds = hl.vds.split_multi(vds)
        mt = split_vds.variant_data

        weights_ht_processed = weights_table.annotate(
            alleles=hl.if_else(
                ref_is_effect_allele,
                [weights_table.effect_allele, weights_table.noneffect_allele],
                [weights_table.noneffect_allele, weights_table.effect_allele],
            ),
            weight=hl.if_else(
                ref_is_effect_allele,
                -weights_table.weight,
                weights_table.weight,
            ),
        ).key_by('locus', 'alleles')

        mt = mt.annotate_rows(
            weights_info=weights_ht_processed[mt.row_key]
        )
        mt = mt.filter_rows(hl.is_defined(mt.weights_info))

    with _log_timing(
        "Planning: Calculating per-variant dosage",
        detailed_timings,
    ):
        # After splitting, LGT is converted to GT, so we can
        # directly and safely use the built-in dosage calculator.
        # See the source code for `hl.vds.split_multi` for details.
        mt = mt.annotate_entries(dosage=mt.GT.n_alt_alleles())

        return mt


def _prepare_mt_non_split(
    vds: hl.vds.VariantDataset,
    weights_table: hl.Table,
    strict_allele_match: bool,
    detailed_timings: bool,
) -> hl.MatrixTable:
    """
    Prepares a MatrixTable for the non-split PRS calculation path.

    This function takes an interval-filtered VDS and joins it with the
    weights table using a locus-based key. It optionally performs a strict
    allele match to handle allele orientation and then calculates dosage
    using the custom multi-allelic dosage function.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        The interval-filtered Variant Dataset.
    weights_table : hail.Table
        The chunk of the weights table.
    strict_allele_match : bool
        If True, performs a robust check for allele correspondence.
    detailed_timings : bool
        If True, logs the duration of computational steps.

    Returns
    -------
    hail.MatrixTable
        A prepared MatrixTable, filtered and annotated with `weights_info`
        and `dosage` for the specified effect allele.
    """
    mt = vds.variant_data

    with _log_timing(
        "Planning: Annotating variants with weights", detailed_timings
    ):
        mt = mt.annotate_rows(weights_info=weights_table[mt.locus])
        mt = mt.filter_rows(hl.is_defined(mt.weights_info))

    if strict_allele_match:
        with _log_timing("Performing strict allele match", detailed_timings):
            alt_alleles = hl.set(mt.alleles[1:])
            ref_allele = mt.alleles[0]
            effect = mt.weights_info.effect_allele
            noneffect = mt.weights_info.noneffect_allele
            is_valid_pair = (
                (
                    (effect == ref_allele)
                    & alt_alleles.contains(noneffect)
                )
                | (
                    (noneffect == ref_allele)
                    & alt_alleles.contains(effect)
                )
            )
            mt = mt.filter_rows(is_valid_pair)

    with _log_timing(
        "Planning: Calculating per-variant dosage",
        detailed_timings,
    ):
        mt = mt.annotate_entries(dosage=_calculate_dosage(mt))

    return mt


def _calculate_prs_chunk(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    include_n_shared_loci: bool,
    sample_id_col: str,
    split_multi: bool,
    ref_is_effect_allele: bool,
    strict_allele_match: bool,
    detailed_timings: bool,
) -> hl.Table:
    """
    Calculates a Polygenic Risk Score (PRS) for a single chunk of variants.

    This function serves as the core computational engine for a single chunk of
    a larger weights table, with an option to handle multi-allelic variants
    by splitting them into bi-allelic representations.

    Parameters
    ----------
    weights_table : hail.Table
        A pre-filtered chunk of the main weights table, keyed by 'locus'.
    vds : hail.vds.VariantDataset
        The Variant Dataset containing the genetic data for the cohort.
    include_n_shared_loci : bool
        If True, the output table will include a column with the total number
        of loci used in the score calculation.
    sample_id_col : str
        The desired name for the sample ID column in the final output table.
    split_multi : bool, default False
        If True, splits multi-allelic variants in the VDS into bi-allelic
        variants before calculation.
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
    detailed_timings : bool
        If True, logs the duration of key computational steps.

    Returns
    -------
    hail.Table
        A Hail Table containing the partial PRS results for the chunk.
    """
    with _log_timing(
        "Planning: Filtering VDS to variants in weights table chunk",
        detailed_timings,
    ):
        intervals_ht = (
            weights_table.select(
                locus_interval=hl.interval(
                    weights_table.locus,
                    weights_table.locus,
                    includes_end=True,
                )
            )
            .key_by('locus_interval')
            .distinct()
        )
        vds = hl.vds.filter_intervals(vds, intervals_ht, keep=True)

    if split_multi:
        mt = _prepare_mt_split(
            vds=vds,
            weights_table=weights_table,
            ref_is_effect_allele=ref_is_effect_allele,
            detailed_timings=detailed_timings,
        )
    else:
        mt = _prepare_mt_non_split(
            vds=vds,
            weights_table=weights_table,
            strict_allele_match=strict_allele_match,
            detailed_timings=detailed_timings,
        )

    # Chunks aggregation
    prs_table = mt.select_cols(
        prs=hl.agg.sum(mt.dosage * mt.weights_info.weight)
    ).cols()

    if include_n_shared_loci:
        with _log_timing("Computing shared loci count", detailed_timings):
            # Using hl.agg.count() within the `select_cols` block won't work
            # since homozygous reference are set to missing while `agg.count`
            # counts the number of rows for which that specific sample has a
            # non-missing genotype calls.
            # This is two-pass approach and thus less performant.
            n_shared_loci = mt.count_rows()
            logger.info("%d loci in common in this chunk.", n_shared_loci)
            prs_table = prs_table.annotate(n_shared_loci=n_shared_loci)

    prs_table = prs_table.rename({'s': sample_id_col})
    # Dropping all global annotations by passing no arguments to select_globals().
    # This is intentional to reduce memory usage.
    return prs_table.select_globals()


def calculate_prs(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    output_path: str,
    chunk_size: int = 20000,
    samples_to_keep: typing.Optional[typing.Any] = None,
    weight_col_name: str = 'weight',
    log_transform_weight: bool = False,
    include_n_shared_loci: bool = False,
    sample_id_col: str = 'person_id',
    split_multi: bool = True,
    ref_is_effect_allele: bool = False,
    strict_allele_match: bool = True,
    detailed_timings: bool = False,
) -> typing.Optional[str]:
    """
    Calculates PRS by chunking weights, checkpointing each chunk, aggregating
    the results, and exporting the final table to a text file.

    This function is optimized for the All of Us VDS, where multi-allelic
    variants are not pre-split. It is compatible with different VDS versions by
    handling both global ('GT') and local ('LGT'/'LA') genotype encoding
    schemes.

    Notes
    -----
    By default (split_multi=True), this function prioritizes robustness over
    performance by splitting multi-allelic variants. This is the best method
    for handling the complex variants common in WGS data.

    This split_multi process includes creating a minimal representation for
    variants. For example, for a variant chr1:10075251 A/G in the weights
    table, split_multi can intelligently match it to a complex indel in the VDS
    (e.g., alleles=['AGGGC', 'A', 'GGGGC']) by simplifying the VDS
    representation to its minimal form (['A', 'G']) for 'AGGGC' -> 'GGGGC'.

    The non-split path (`split_multi=False`) is a faster but less robust
    alternative. It relies on a direct string comparison of alleles and will
    fail to match the complex variant described above. Furthermore, if the
    weights table contains multiple entries for the same locus, the non-split
    path will arbitrarily select only one of them. This "power-user" option
    should only be used if you are certain that both your VDS and weights table
    contain only simple, well-matched, bi-allelic variants.

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
        The Variant Dataset containing the genetic data.
    output_path : str
        A GCS path (starting with 'gs://') to write the final output file.
    chunk_size : int, default 20000
        The number of variants to include in each processing chunk.
    samples_to_keep : hail.Table, list, set, tuple, int, or str, optional
        A collection of sample IDs to keep. If provided, PRS will only be
        calculated for these samples. Numeric IDs will be converted to strings.
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
        If True, splits multi-allelic variants in the VDS into bi-allelic
        variants before calculation.
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

    Returns
    -------
    str or None
        The path to the final output file. Returns None if no results are
        generated. The output file is a tab-separated text file with the
        following columns:
        - A sample identifier column (named according to `sample_id_col`).
        - 'prs': The calculated Polygenic Risk Score.
        - 'n_shared_loci' (optional): The number of variants used to
          calculate the score, included if `include_n_shared_loci` is True.

    Raises
    ------
    ValueError
        If `output_path` or `temp_path` are not valid GCS paths, or if the
        `weights_table` is empty after validation.

    """
    timer = SimpleTimer()
    with timer:
        if not output_path.startswith('gs://'):
            raise ValueError(
                "The 'output_path' must be a Google Cloud Storage (GCS) "
                "path, starting with 'gs://'."
            )

        logger.info(
            "Starting PRS calculation. Final result will be at: %s",
            output_path,
        )

        if samples_to_keep is not None:
            with _log_timing(
                "Filtering to specified samples", detailed_timings
            ):
                samples_ht = _prepare_samples_to_keep(samples_to_keep)
                vds = hl.vds.filter_samples(vds, samples_ht)

        with _log_timing(
            "Preparing and chunking full weights table", detailed_timings
        ):
            full_weights_table = _validate_and_prepare_weights_table(
                weights_table, weight_col_name, log_transform_weight
            )
            total_variants = full_weights_table.count()
            if total_variants == 0:
                raise ValueError(
                    "Weights table is empty after validation."
                )

            n_chunks = ceil(total_variants / chunk_size)
            logger.info(
                "Total variants: %d, Number of chunks: %d",
                total_variants,
                n_chunks,
            )

            full_weights_table = full_weights_table.add_index()
            full_weights_table = full_weights_table.annotate(
                chunk_id=hl.int(full_weights_table.idx / chunk_size)
            )
            full_weights_table = full_weights_table.key_by('locus')

        partial_dfs = []
        for i in range(n_chunks):
            with _log_timing(
                f"Processing chunk {i + 1}/{n_chunks}",
                True,
            ):
                # Use .persist() to avoid recomputation of the same chunk in
                # _calculate_prs_chunk, specifically during:
                # 1. Creation of interval_ht
                # 2. Annotating rows with PRS weight information
                weights_chunk = full_weights_table.filter(
                    full_weights_table.chunk_id == i
                ).persist()

                chunk_prs_table = _calculate_prs_chunk(
                    weights_table=weights_chunk,
                    vds=vds,
                    include_n_shared_loci=include_n_shared_loci,
                    sample_id_col=sample_id_col,
                    split_multi=split_multi,
                    ref_is_effect_allele=ref_is_effect_allele,
                    strict_allele_match=strict_allele_match,
                    detailed_timings=detailed_timings,
                )

                partial_dfs.append(chunk_prs_table.to_pandas())

        if not partial_dfs:
            logger.warning("No PRS results were generated.")
            return None

        # Aggregate the results from all chunks in master node
        with _log_timing(
            f"Aggregating chunks and exporting to {output_path}",
            True,
        ):
            combined_df = pd.concat(partial_dfs, ignore_index=True)
            final_df = combined_df.groupby(sample_id_col).sum()

        # Export the final pandas data frame to a csv file
        with _log_timing(
                f"Exporting final result to {output_path}", detailed_timings
        ):
            with hfs.open(output_path, 'w') as f:
                final_df.to_csv(f, sep='\t', index=True, header=True)

    # Report the total time using the duration captured by the context manager
    logger.info(
        "PRS calculation complete. Total time: %.2f seconds.", timer.duration
    )
    return output_path


# This function is a legacy version of the PRS calculation that merges chunks
# using Hail Table joins.
def calculate_prs2(
    weights_table: hl.Table,
    vds: hl.vds.VariantDataset,
    output_path: str,
    temp_path: typing.Optional[str] = None,
    chunk_size: int = 20000,
    samples_to_keep: typing.Optional[typing.Any] = None,
    weight_col_name: str = 'weight',
    log_transform_weight: bool = False,
    include_n_shared_loci: bool = False,
    sample_id_col: str = 'person_id',
    split_multi: bool = True,
    ref_is_effect_allele: bool = False,
    strict_allele_match: bool = True,
    detailed_timings: bool = False,
) -> typing.Optional[str]:
    """
    This function is a legacy version of `calculate_prs` that merges chunks
    using Hail Table joins, with temporary checkpointing.

    Calculates PRS by chunking weights, checkpointing each chunk, aggregating
    the results, and exporting the final table to a text file.

    This function is optimized for the All of Us VDS, where multi-allelic
    variants are not pre-split. It is compatible with different VDS versions by
    handling both global ('GT') and local ('LGT'/'LA') genotype encoding
    schemes.

    Notes
    -----
    By default (split_multi=True), this function prioritizes robustness over
    performance by splitting multi-allelic variants. This is the best method
    for handling the complex variants common in WGS data.

    This split_multi process includes creating a minimal representation for
    variants. For example, for a variant chr1:10075251 A/G in the weights
    table, split_multi can intelligently match it to a complex indel in the VDS
    (e.g., alleles=['AGGGC', 'A', 'GGGGC']) by simplifying the VDS
    representation to its minimal form (['A', 'G']) for 'AGGGC' -> 'GGGGC'.

    The non-split path (`split_multi=False`) is a faster but less robust
    alternative. It relies on a direct string comparison of alleles and will
    fail to match the complex variant described above. Furthermore, if the
    weights table contains multiple entries for the same locus, the non-split
    path will arbitrarily select only one of them. This "power-user" option
    should only be used if you are certain that both your VDS and weights table
    contain only simple, well-matched, bi-allelic variants.

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
        The Variant Dataset containing the genetic data.
    output_path : str
        A GCS path (starting with 'gs://') to write the final output file.
    temp_path : str, optional
        A GCS path for storing temporary checkpoint files. If not provided,
        it will default to the 'WORKSPACE_BUCKET' environment variable.
    chunk_size : int, default 20000
        The number of variants to include in each processing chunk.
    samples_to_keep : hail.Table, list, set, tuple, int, or str, optional
        A collection of sample IDs to keep. If provided, PRS will only be
        calculated for these samples. Numeric IDs will be converted to strings.
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
        If True, splits multi-allelic variants in the VDS into bi-allelic
        variants before calculation.
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

    Returns
    -------
    str or None
        The path to the final output file. Returns None if no results are
        generated. The output file is a tab-separated text file with the
        following columns:
        - A sample identifier column (named according to `sample_id_col`).
        - 'prs': The calculated Polygenic Risk Score.
        - 'n_shared_loci' (optional): The number of variants used to
          calculate the score, included if `include_n_shared_loci` is True.

    Raises
    ------
    ValueError
        If `output_path` or `temp_path` are not valid GCS paths, or if the
        `weights_table` is empty after validation.
    """
    timer = SimpleTimer()
    with timer:
        if not output_path.startswith('gs://'):
            raise ValueError(
                "The 'output_path' must be a Google Cloud Storage (GCS) "
                "path, starting with 'gs://'."
            )

        temp_path = temp_path or os.getenv('WORKSPACE_BUCKET')
        if not temp_path or not temp_path.startswith('gs://'):
            raise ValueError(
                "A valid GCS path (starting with 'gs://') for temporary files "
                "must be provided via the 'temp_path' argument or the "
                "'WORKSPACE_BUCKET' environment variable."
            )

        logger.info(
            "Starting PRS calculation. Final result will be at: %s",
            output_path,
        )

        temp_dir = f"{temp_path.rstrip('/')}/prs_temp_{uuid4()}"
        logger.info("Using temporary checkpointing directory: %s", temp_dir)

        try:
            if samples_to_keep is not None:
                with _log_timing(
                    "Filtering to specified samples", detailed_timings
                ):
                    samples_ht = _prepare_samples_to_keep(samples_to_keep)
                    vds = hl.vds.filter_samples(vds, samples_ht)

            with _log_timing(
                "Preparing and chunking full weights table", detailed_timings
            ):
                full_weights_table = _validate_and_prepare_weights_table(
                    weights_table, weight_col_name, log_transform_weight
                )
                total_variants = full_weights_table.count()
                if total_variants == 0:
                    raise ValueError(
                        "Weights table is empty after validation."
                    )

                n_chunks = ceil(total_variants / chunk_size)
                logger.info(
                    "Total variants: %d, Number of chunks: %d",
                    total_variants,
                    n_chunks,
                )

                full_weights_table = full_weights_table.add_index()
                full_weights_table = full_weights_table.annotate(
                    chunk_id=hl.int(full_weights_table.idx / chunk_size)
                )
                full_weights_table = full_weights_table.key_by('locus')

            partial_results = []
            for i in range(n_chunks):
                with _log_timing(
                    f"Processing and checkpointing chunk {i + 1}/{n_chunks}",
                    True,
                ):
                    # Use .persist() to avoid recomputation of the same chunk in
                    # _calculate_prs_chunk, specifically during:
                    # 1. Creation of interval_ht
                    # 2. Annotating rows with PRS weight information
                    weights_chunk = full_weights_table.filter(
                        full_weights_table.chunk_id == i
                    ).persist()
                    checkpoint_path = f"{temp_dir}/chunk_{i}.ht"

                    chunk_prs_table = _calculate_prs_chunk(
                        weights_table=weights_chunk,
                        vds=vds,
                        include_n_shared_loci=include_n_shared_loci,
                        sample_id_col=sample_id_col,
                        split_multi=split_multi,
                        ref_is_effect_allele=ref_is_effect_allele,
                        strict_allele_match=strict_allele_match,
                        detailed_timings=detailed_timings,
                    )

                    processed_chunk = chunk_prs_table.checkpoint(
                        checkpoint_path, overwrite=True
                    )
                    partial_results.append(processed_chunk)

            if not partial_results:
                logger.warning(
                    "No PRS results were generated. "
                    "No output file will be created."
                )
                return None

            with _log_timing(
                f"Aggregating chunks and exporting to {output_path}",
                True,
            ):
                final_prs_table = partial_results[0].key_by(sample_id_col)
                for i in range(1, len(partial_results)):
                    next_chunk_table = partial_results[i].key_by(
                        sample_id_col
                    )
                    final_prs_table = final_prs_table.join(
                        next_chunk_table, how='outer'
                    )

                    current_prs = hl.coalesce(final_prs_table.prs, 0)
                    next_prs = hl.coalesce(final_prs_table.prs_1, 0)
                    final_prs_table = final_prs_table.annotate(
                        prs=current_prs + next_prs
                    )

                    if include_n_shared_loci:
                        current_loci = hl.coalesce(
                            final_prs_table.n_shared_loci, 0
                        )
                        next_loci = hl.coalesce(
                            final_prs_table.n_shared_loci_1, 0
                        )
                        final_prs_table = final_prs_table.annotate(
                            n_shared_loci=current_loci + next_loci
                        ).key_by(sample_id_col)

                    cols_to_keep = ['prs']
                    if include_n_shared_loci:
                        cols_to_keep.append('n_shared_loci')
                    final_prs_table = final_prs_table.select(*cols_to_keep)

                final_prs_table.export(
                    output_path, header=True, delimiter='\t'
                )

        finally:
            logger.info("Cleaning up temporary directory: %s", temp_dir)
            if hfs.is_dir(temp_dir):
                hfs.rmtree(temp_dir)

    # Report the total time using the duration captured by the context manager
    logger.info(
        "PRS calculation complete. Total time: %.2f seconds.", timer.duration
    )
    return output_path


# def calculate_prs_batch(
#     weights_map: typing.Dict[str, hl.Table],
#     vds: hl.vds.VariantDataset,
#     samples_to_keep: typing.Optional[typing.Any] = None,
#     weight_col_name: str = 'weight',
#     log_transform_weight: bool = False,
#     strict_allele_match: bool = True,
#     include_n_shared_loci: bool = False,
#     sample_id_col: str = 'person_id',
#     detailed_timings: bool = True
# ) -> hl.Table:
#     """
#     Calculates multiple Polygenic Risk Scores (PRS) in a single pass.

#     This function is highly efficient for calculating several PRS on the same
#     cohort, as it filters the VDS only once for all combined variants.

#     Note: it is optimized for performance on the All of Us VDS, where
#     multi-allelic variants are not pre-split. It computes allele dosages
#     directly from the multi-allelic VDS, avoiding the computationally expensive
#     `split_multi` step. It is compatible with different VDS versions by
#     handling both global ('GT') and local ('LGT'/'LA') genotype encoding
#     schemes.

#     To achieve this, the function joins the VDS and the weights table using
#     only the genomic `locus` as a key. A critical consequence of this design is
#     that it assumes there is only **one** variant entry per locus in the
#     `weights_table`. If the weights table contains multiple entries for the
#     same locus, Hail will use the first matching entry it encounters, which can
#     lead to an incorrect PRS calculation. If you prefer additional robustness
#     over speed, consider using the `calculate_prs_split_batch` function.

#     Parameters
#     ----------
#     weights_map : dict[str, hail.Table]
#         A dictionary where keys are the desired PRS names (e.g., 'CAD_prs')
#         and values are the corresponding Hail Tables of weights. Each table
#         must contain the required columns ('chr', 'pos', 'effect_allele',
#         'noneffect_allele', and a weight column).
#     vds : hail.vds.VariantDataset
#         A Variant Dataset object containing the genetic data.
#     samples_to_keep : hail.Table, list, set, tuple, int, or str, optional
#         A collection of sample IDs to keep.
#     weight_col_name : str, default 'weight'
#         The name of the column in `weights_table` that contains the effect
#         weights.
#     log_transform_weight : bool, default False
#         If True, applies a natural log transformation to the weight column.
#         This should be used when weights are provided as odds ratios (OR).
#     strict_allele_match : bool, default True
#         If True, validates that one of the alleles from weight table is the
#         reference and the other is a valid alternate allele in the VDS.
#     include_n_shared_loci : bool, default False
#         If True, includes a column for each score with the total number of
#         shared loci used. Note that this option can impact performance.
#     sample_id_col : str, default 'person_id'
#         The desired name for the sample ID column in the final output table.
#     detailed_timings : bool, default True
#         If True, logs the duration of each major computational step.

#     Returns
#     -------
#     hail.Table
#         A "wide" Hail Table with one column per calculated PRS and one
#         column for the number of shared loci for each score.

#     """
#     logger.info(
#         f"Starting batch PRS calculation for {len(weights_map)} scores..."
#     )

#     if samples_to_keep is not None:
#         with _log_timing("Filtering to specified samples", detailed_timings):
#             samples_ht = _prepare_samples_to_keep(samples_to_keep)
#             vds = hl.vds.filter_samples(vds, samples_ht)

#     prepared_weights = {}
#     all_loci = []
#     with _log_timing("Preparing all weights tables", detailed_timings):
#         for score_name, weights_table in weights_map.items():
#             prepared_table = _validate_and_prepare_weights_table(
#                 weights_table, weight_col_name, log_transform_weight
#             )
#             prepared_weights[score_name] = prepared_table
#             all_loci.append(prepared_table.key_by('locus').select())

#     with _log_timing("Filtering VDS to all variants", detailed_timings):
#         loci_to_keep = hl.Table.union(*all_loci).key_by('locus')
#         count = loci_to_keep.count()
#         logger.info(f"Found {count} total unique loci across all scores.")
#         mt = vds.variant_data.semi_join_rows(loci_to_keep)

#     with _log_timing("Annotating variants with all weights", detailed_timings):
#         annotation_exprs = {
#             f'weights_info_{score_name}': prepared_weights[score_name][mt.locus]
#             for score_name in weights_map
#         }
#         mt = mt.annotate_rows(**annotation_exprs)

#     if include_n_shared_loci:
#         logger.info("Calculating scores and n_shared_loci.")
#         if strict_allele_match:
#             with _log_timing("Flagging valid loci", detailed_timings):
#                 vds_alleles = hl.set(mt.alleles)
#                 ref_allele = mt.alleles[0]
#                 validity_annotations = {}
#                 for score_name in weights_map:
#                     weights_info = mt[f'weights_info_{score_name}']
#                     is_present = hl.is_defined(weights_info)
#                     effect = weights_info.effect_allele
#                     noneffect = weights_info.noneffect_allele
#                     is_valid_pair = (
#                         (
#                             (effect == ref_allele)
#                             & vds_alleles.contains(noneffect)
#                         ) | (
#                             (noneffect == ref_allele)
#                             & vds_alleles.contains(effect)
#                         )
#                     ) & (effect != noneffect)
#                     validity_annotations[f'{score_name}_is_valid'] = hl.if_else(
#                         is_present, is_valid_pair, False
#                     )
#                 mt = mt.annotate_rows(**validity_annotations)

#         with _log_timing("Pre-calculating n_shared_loci", detailed_timings):
#             counts_by_score = mt.aggregate_rows(
#                 hl.struct(**{
#                     f'{score_name}_n_loci': hl.agg.count_where(
#                         mt[f'{score_name}_is_valid']
#                     )
#                     for score_name in weights_map
#                 })
#             )

#         with _log_timing("Calculating PRS scores", detailed_timings):
#             aggregators = {}
#             for score_name in weights_map:
#                 weights_info = mt[f'weights_info_{score_name}']
#                 variant_filter = mt[f'{score_name}_is_valid']
#                 dosage_expr = _calculate_dosage(mt, score_name)
#                 aggregators[f'{score_name}_prs'] = hl.agg.filter(
#                     variant_filter,
#                     hl.agg.sum(dosage_expr * weights_info.weight)
#                 )
#             mt = mt.annotate_cols(**aggregators)

#         logger.info("Finalizing batch PRS table...")
#         prs_table = mt.cols()
#         final_select_exprs = {}
#         for score_name in weights_map.keys():
#             final_select_exprs[score_name] = prs_table[f'{score_name}_prs']
#             final_select_exprs[f'{score_name}_n_loci'] = (
#                 counts_by_score[f'{score_name}_n_loci']
#             )
#         prs_table = prs_table.select(**final_select_exprs)

#     else:
#         logger.info("Calculating scores only (n_shared_loci not requested).")
#         if strict_allele_match:
#             with _log_timing("Performing strict allele match", detailed_timings):
#                 vds_alleles = hl.set(mt.alleles)
#                 ref_allele = mt.alleles[0]
#                 new_annotations = {}
#                 for score_name in weights_map:
#                     weights_info = mt[f'weights_info_{score_name}']
#                     effect = weights_info.effect_allele
#                     noneffect = weights_info.noneffect_allele
#                     is_valid_pair = (
#                         (
#                             (effect == ref_allele)
#                             & vds_alleles.contains(noneffect)
#                         ) | (
#                             (noneffect == ref_allele)
#                             & vds_alleles.contains(effect)
#                         )
#                     ) & (effect != noneffect)
#                     new_annotations[f'weights_info_{score_name}'] = hl.if_else(
#                         is_valid_pair,
#                         weights_info,
#                         hl.missing(weights_info.dtype)
#                     )
#                 mt = mt.annotate_rows(**new_annotations)

#         with _log_timing("Calculating PRS scores", detailed_timings):
#             aggregators = {}
#             for score_name in weights_map:
#                 weights_info = mt[f'weights_info_{score_name}']
#                 variant_filter = hl.is_defined(weights_info)
#                 dosage_expr = _calculate_dosage(mt, score_name)
#                 aggregators[f'{score_name}_prs'] = hl.agg.filter(
#                     variant_filter,
#                     hl.agg.sum(dosage_expr * weights_info.weight)
#                 )
#             mt = mt.annotate_cols(**aggregators)

#         logger.info("Finalizing batch PRS table...")
#         prs_table = mt.cols()
#         prs_table = prs_table.select(
#             *[f'{s}_prs' for s in weights_map.keys()]
#         )
#         prs_table = prs_table.rename(
#             {f'{s}_prs': s for s in weights_map.keys()}
#         )

#     prs_table = prs_table.rename({'s': sample_id_col})
#     logger.info("Batch PRS calculation complete.")
#     return prs_table
