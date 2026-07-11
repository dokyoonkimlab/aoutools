"""Internal utilities for the PRS calculator"""

# See the note in _config.py: mocked hail (docs builds) cannot evaluate the
# PEP 604 `|` in these annotations, so defer them.
from __future__ import annotations

import logging
from math import ceil

import hail as hl

from ._config import PRSConfig
from ._utils import (
    _log_timing,
    _standardize_chromosome_column,
)

logger = logging.getLogger(__name__)


def _prepare_samples_to_keep(
    samples: hl.Table | list | set | tuple | int | float | str,
) -> hl.Table:
    """
    Converts a flexible list of samples into a keyed Hail Table.

    Accepts various common Python collection types (list, set, tuple) or single
    values (int, float, str) and converts them into a standardized Hail Table
    keyed by the string column 's'. This table can then be used to filter
    Hail datasets by sample ID.

    Parameters
    ----------
    samples : hail.Table, list, set, tuple, int, float, or str
        A collection or single value of sample IDs to prepare. Numeric types
        will be converted to strings for consistent keying.

    Returns
    -------
    hail.Table
        A Hail Table keyed by 's' containing sample IDs as strings.

    Raises
    ------
    TypeError
        If the input `samples` object is not one of the supported types.
    """
    if isinstance(samples, hl.Table):
        return samples

    if isinstance(samples, (int, float, str)):
        sample_list = [str(samples)]
    elif isinstance(samples, (list, set, tuple)):
        sample_list = [str(s) for s in samples]
    else:
        raise TypeError(
            f"Unsupported type for samples_to_keep: {type(samples)}."
        )

    samples_ht = hl.Table.parallelize(
        [{"s": s} for s in sample_list], hl.tstruct(s=hl.tstr)
    )
    return samples_ht.key_by("s")


def _validate_and_prepare_weights_table(
    weights_table: hl.Table, config: PRSConfig
) -> hl.Table:
    """
    Validates and prepares a single weights table for PRS calculation.

    Ensures the table has the required columns with correct types, standardizes
    the chromosome format, handles the weight column (renaming and optional
    log transformation), and keys the table by locus for joining with the VDS.

    Parameters
    ----------
    weights_table : hail.Table
        An input weights table. Must contain columns 'chr', 'pos',
        'effect_allele', 'noneffect_allele', and a weight column.
    config : PRSConfig
        A configuration object specifying parameters such as `weight_col_name`
        and `log_transform_weight`.

    Returns
    -------
    hail.Table
        A validated, standardized, and keyed Hail Table.

    Raises
    ------
    TypeError
        If `weight_col_name` is missing, or any required column is missing, or
        has incorrect data type.

    See also
    --------
    PRSConfig : A configuration class that holds parameters for PRS calculation.
    """
    if config.weight_col_name not in weights_table.row:
        raise TypeError(
            f"Specified weight column '{config.weight_col_name}' not found "
            f"in table."
        )
    weights_table = weights_table.rename({config.weight_col_name: "weight"})

    required_cols = {
        "chr": hl.tstr,
        "pos": hl.tint32,
        "effect_allele": hl.tstr,
        "noneffect_allele": hl.tstr,
        "weight": hl.tfloat64,
    }
    for col, expected_type in required_cols.items():
        if col not in weights_table.row:
            raise TypeError(
                f"Weights table is missing required column: '{col}'."
            )
        if weights_table[col].dtype != expected_type:
            raise TypeError(f"Column '{col}' has incorrect type.")

    weights_table = _standardize_chromosome_column(weights_table)

    if config.log_transform_weight:
        weights_table = weights_table.annotate(
            weight=hl.log(weights_table.weight)
        )

    weights_table = weights_table.annotate(
        locus=hl.locus(
            weights_table.chr, weights_table.pos, reference_genome="GRCh38"
        )
    )
    weights_table = weights_table.key_by("locus")
    return weights_table.select("effect_allele", "noneffect_allele", "weight")


def _orient_weights_for_split(ht: hl.Table, config: PRSConfig) -> hl.Table:
    """
    Orients alleles and weights for a split-multi join.

    Constructs a canonical `[ref, alt]` allele representation for the join key
    and adjusts the weights so that they always correspond to the alternate
    allele. This is important to ensure consistency in allele orientation for
    PRS calculation.

    Parameters
    ----------
    ht : hail.Table
        A Hail table containing 'effect_allele', 'noneffect_allele', 'weight',
        and 'locus' fields.
    config : PRSConfig
        A configuration object specifying if reference allele is effect allele
        (`ref_is_effect_allele`).

    Returns
    -------
    hail.Table
        A table keyed by 'locus' and 'alleles', with adjusted 'weight' to align
        with the alternate allele.

    See also
    --------
    PRSConfig : A configuration class that specifies allele orientation
        settings.
    """
    return ht.annotate(
        alleles=hl.if_else(
            config.ref_is_effect_allele,
            [ht.effect_allele, ht.noneffect_allele],
            [ht.noneffect_allele, ht.effect_allele],
        ),
        weight=hl.if_else(config.ref_is_effect_allele, -ht.weight, ht.weight),
    ).key_by("locus", "alleles")


def _prepare_weights_for_chunking(
    weights_table: hl.Table, config: PRSConfig, validate_table: bool = True
) -> tuple[hl.Table, int]:
    """
    Prepares and annotates a weights table for chunked processing.

    This helper function takes a raw weights table, optionally validates it,
    and assigns each row a `chunk_id` based on the configured chunk size.
    It enables iterative PRS calculation by partitioning the input into
    manageable chunks.

    Parameters
    ----------
    weights_table : hail.Table
        Raw input weights table from the user.
    weight_col_name : str
        Column name containing effect weights.
    config : PRSConfig
        Configuration object with `chunk_size`, `log_transform_weight`, and
        `detailed_timings` settings.
    validate_table : bool, default=True
        If True, validates and preprocesses the weights table using
        `_validate_and_prepare_weights_table`. If False, assumes the table is
        already validated. This is useful for batch processing where validation
        occurs upstream.

    Returns
    -------
    tuple[hail.Table, int]
        Tuple containing:
        - A validated and annotated weights table with a `chunk_id` column.
        - The number of chunks the table was divided into.

    Raises
    ------
    ValueError
        If the input table is empty after validation.
    """
    with _log_timing(
        "Preparing and analyzing weights table", config.detailed_timings
    ):
        if validate_table:
            weights_table = _validate_and_prepare_weights_table(
                weights_table=weights_table, config=config
            )

        total_variants = weights_table.count()
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
        weights_table = weights_table.add_index()
        weights_table = weights_table.annotate(
            chunk_id=hl.int(weights_table.idx / effective_chunk_size)
        )
        # Note: add_index() preserves existing keys (e.g., locus),
        # so no need to re-key explicitly

        return weights_table, n_chunks


def _create_1bp_intervals(table_chunk: hl.Table) -> hl.Table:
    """
    Creates a table of 1-base-pair (1-bp) intervals from a Hail Table keyed by
    locus.

    Useful for interval-based joins or filtering operations, where each variant
    is represented as a genomic interval spanning exactly one position
    (i.e., [locus, locus], inclusive).

    Parameters
    ----------
    table_chunk : hail.Table
        A Hail Table containing a 'locus' field of type `locus<Locus>`.

    Returns
    -------
    hail.Table
        A new Table keyed by 1-bp interval around each locus.
    """
    return table_chunk.select(
        interval=hl.interval(
            table_chunk.locus, table_chunk.locus, includes_end=True
        )
    ).key_by("interval")
