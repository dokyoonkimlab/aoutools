"""Internal utilities for the PRS calculator"""

import typing
import logging
from math import ceil
import hail as hl
from ._utils import (
    _log_timing,
    _standardize_chromosome_column,
)
from ._config import PRSConfig

logger = logging.getLogger(__name__)


def _prepare_samples_to_keep(
    samples: typing.Union[hl.Table, list, set, tuple, int, str]
) -> hl.Table:
    """
    Converts a flexible list of samples into a keyed Hail Table.

    This helper function provides flexibility by accepting various common
    Python collection types (list, set, tuple) or single values (int, str)
    and converting them into a standardized Hail Table with a string key 's',
    which is required for filtering Hail objects.

    Parameters
    ----------
    samples : hail.Table, list, set, tuple, int, or str
        The collection of sample IDs to prepare.

    Returns
    -------
    hail.Table
        A Hail Table keyed by 's' containing the sample IDs as strings.

    Raises
    ------
    TypeError
        If the input `samples` object is not one of the supported types.
    """
    if isinstance(samples, hl.Table):
        return samples

    sample_list = []
    if isinstance(samples, (int, float, str)):
        sample_list = [str(samples)]
    elif isinstance(samples, (list, set, tuple)):
        sample_list = [str(s) for s in samples]
    else:
        raise TypeError(f"Unsupported type for samples_to_keep: {type(samples)}.")

    samples_ht = hl.Table.parallelize(
        [{'s': s} for s in sample_list], hl.tstruct(s=hl.tstr)
    )
    return samples_ht.key_by('s')


def _validate_and_prepare_weights_table(
    weights_table: hl.Table,
    config: PRSConfig
) -> hl.Table:
    """
    Validates and prepares a single weights table for PRS calculation.

    This function ensures the table has the required columns with the correct
    types, standardizes the chromosome format, handles the weight column,
    and keys the table by locus for joining with the VDS.

    Parameters
    ----------
    weights_table : hail.Table
        The input weights table. Must contain 'chr', 'pos', 'effect_allele',
        'noneffect_allele', and a weight column.
    config : PRSConfig
        A configuration object containing settings like `weight_col_name` and
        `log_transform_weight`.

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
    if config.weight_col_name not in weights_table.row:
        raise TypeError(
            f"Specified weight column '{config.weight_col_name}' not found "
            f"in table."
        )
    weights_table = weights_table.rename({config.weight_col_name: 'weight'})

    required_cols = {
        'chr': hl.tstr,
        'pos': hl.tint32,
        'effect_allele': hl.tstr,
        'noneffect_allele': hl.tstr,
        'weight': hl.tfloat64,
    }
    for col, expected_type in required_cols.items():
        if col not in weights_table.row:
            raise TypeError(f"Weights table is missing required column: '{col}'.")
        if weights_table[col].dtype != expected_type:
            raise TypeError(f"Column '{col}' has incorrect type.")

    weights_table = _standardize_chromosome_column(weights_table)

    if config.log_transform_weight:
        weights_table = weights_table.annotate(
            weight=hl.log(weights_table.weight)
        )

    weights_table = weights_table.annotate(
        locus=hl.locus(
            weights_table.chr,
            weights_table.pos,
            reference_genome='GRCh38'
        )
    )
    weights_table = weights_table.key_by('locus')
    return weights_table.select('effect_allele', 'noneffect_allele', 'weight')


def _orient_weights_for_split(
    ht: hl.Table,
    config: PRSConfig
) -> hl.Table:
    """
    Orients alleles and weight for a split-multi join.

    Creates a canonical [ref, alt] representation for the join key and
    adjusts the weight to always correspond to the alternate allele.
    """
    return ht.annotate(
        alleles=hl.if_else(
            config.ref_is_effect_allele,
            [ht.effect_allele, ht.noneffect_allele],
            [ht.noneffect_allele, ht.effect_allele]
        ),
        weight=hl.if_else(
            config.ref_is_effect_allele,
            -ht.weight,
            ht.weight
        )
    ).key_by('locus', 'alleles')

def _check_allele_match(
    mt: hl.MatrixTable,
    weights_info: hl.expr.StructExpression
) -> hl.expr.BooleanExpression:
    """
    Returns a boolean expression for a strict allele match.
    This checks if one allele from the weights table is the ref allele in the
    VDS and the other is a valid alt allele.
    """
    alt_alleles = hl.set(mt.alleles[1:])
    ref_allele = mt.alleles[0]
    effect = weights_info.effect_allele
    noneffect = weights_info.noneffect_allele

    is_valid_pair = (
        ((effect == ref_allele) & alt_alleles.contains(noneffect)) |
        ((noneffect == ref_allele) & alt_alleles.contains(effect))
    )

    return is_valid_pair

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


def _prepare_weights_for_chunking(
    weights_table: hl.Table,
    # weight_col_name: str,
    # log_transform_weight: bool,
    # chunk_size: typing.Optional[int],
    # detailed_timings: bool,
    config: PRSConfig,
    validate_table: bool = True
) -> tuple[hl.Table, int]:
    """
    Prepares and annotates a weights table for chunked processing.

    This helper function takes the raw weights table, validates it, calculates
    the number of chunks based on the `chunk_size`, and adds a `chunk_id`
    column to each row. This prepares the table for iterative processing in
    the main PRS calculation loop.

    Parameters
    ----------
    weights_table : hail.Table
        The raw input weights table from the user.
    weight_col_name : str
        The name of the column containing effect weights.
    log_transform_weight : bool
        If True, log-transforms the weight column.
    chunk_size : int, optional
        The desired number of variants per chunk. If None, the entire table
        will be processed as a single chunk.
    detailed_timings : bool
        If True, logs the duration of this preparation step.
    validate_table: bool, default True
        If True, the function calls `_validate_and_prepare_weights_table`. If
        False, this validation step is skipped, assuming the table is already
        prepared. This argument is intended for PRS batch calculations, where
        the table validation has already been performed upstream.

    Returns
    -------
    tuple[hail.Table, int]
        A tuple containing:
        - The fully validated and prepared weights table, now annotated with a
          `chunk_id` for each row.
        - An integer representing the total number of chunks.

    Raises
    ------
    ValueError
        If the `weights_table` is empty after the initial validation and
        filtering steps.
    """
    with _log_timing(
        "Preparing and analyzing weights table", config.detailed_timings
    ):
        if validate_table:
            full_weights_table = _validate_and_prepare_weights_table(
                weights_table=weights_table,
                config=config
            )
        else:
            full_weights_table = weights_table
        total_variants = full_weights_table.count()
        if total_variants == 0:
            raise ValueError("Weights table is empty after validation.")

        effective_chunk_size = config.chunk_size or total_variants
        n_chunks = ceil(total_variants / effective_chunk_size)
        logger.info(
            "Total variants: %d, Number of chunks: %d",
            total_variants,
            n_chunks,
        )

        # Don't use chain (hl.Table.add_index().annotate()) as it is not find
        # the idx at the annotation step due to lazy eval.
        full_weights_table = full_weights_table.add_index()
        full_weights_table = full_weights_table.annotate(
            chunk_id=hl.int(
                full_weights_table.idx / effective_chunk_size
            )
        )
        # Note:add_index does not reset the existing key (locus), so don't have
        # to key_by('locus) again

        return full_weights_table, n_chunks


def _create_1bp_intervals(
        table_chunk: hl.Table
) -> hl.Table:
    """
    Creates a table of 1-bp intervals from a table with a locus key.
    """
    return table_chunk.select(
        interval=hl.interval(
            table_chunk.locus, table_chunk.locus, includes_end=True
        )
    ).key_by('interval')
