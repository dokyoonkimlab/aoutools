import os
import typing
import logging
import hail as hl
from ._utils import (
    _stage_local_file_to_gcs,
    _standardize_chromosome_column,
)

logger = logging.getLogger(__name__)


def _validate_alleles(table: hl.Table) -> hl.Table:
    """
    Filters out rows with invalid alleles (non-ACGT characters).

    This function handles both SNPs and indels by checking that all
    characters in the allele strings are standard DNA bases.

    Parameters
    ----------
    table : hail.Table
        The table to validate.

    Returns
    -------
    hail.Table
        A table with invalid allele rows removed.
    """
    logger.info("Validating allele columns for non-ACGT characters...")
    dna_regex = '^[ACGT]+$'
    initial_count = table.count()
    table = table.filter(
        hl.str(table.effect_allele).matches(dna_regex) &
        hl.str(table.noneffect_allele).matches(dna_regex)
    )
    final_count = table.count()
    n_removed = initial_count - final_count
    if n_removed > 0:
        logger.warning(
            f"Removed {n_removed} variants with invalid alleles (non-ACGT "
            f"characters found)."
        )
    return table


def _check_duplicated_ids(table: hl.Table, file_path: str = "input") -> None:
    """
    Checks for duplicate variants based on genomic identifiers.

    This function creates a temporary variant ID by concatenating the chromosome,
    position, and alleles. It then groups by this ID to identify and
    raise an error if any duplicates are found.

    Parameters
    ----------
    table : hail.Table
        The input table to check. It is expected to have 'chr', 'pos',
        'noneffect_allele', and 'effect_allele' columns.
    file_path : str, optional
        The source file path, used for a more informative error message.
        Defaults to "input".

    Raises
    ------
    ValueError
        If duplicate variants are found in the table.
    """
    logger.info("Checking for duplicate variants based on chr, pos, and alleles...")
    table_with_id = table.annotate(
        _variant_id=hl.str('_').join([
            table.chr, hl.str(table.pos), table.noneffect_allele,
            table.effect_allele
        ])
    )
    id_counts_table = table_with_id.group_by('_variant_id').aggregate(
        n=hl.agg.count()
    )
    duplicate_variants = id_counts_table.filter(id_counts_table.n > 1)
    if duplicate_variants.count() > 0:
        example_duplicates = duplicate_variants.take(5)
        formatted_examples = [d._variant_id for d in example_duplicates]
        raise ValueError(
            f"Duplicate variants found in '{file_path}'. Examples of "
            f"duplicate IDs: {formatted_examples}."
        )
    logger.info("No duplicate variants found.")


def _process_prs_weights_table(
    table: hl.Table, file_path: str, validate_alleles: bool
) -> hl.Table:
    """
    Performs final filtering and validation steps on an imported weights table.

    This function consolidates the shared post-processing logic, including:
    1. Optionally validating that allele columns contain only ACGT characters.
    2. Standardizing the chromosome column to ensure it has a 'chr' prefix.
    3. Filtering out variants with undefined (missing) effect weights.
    4. Checking if the table is empty after all filtering.
    5. Checking for duplicate variants.
    6. Logging the final count of loaded variants.

    Parameters
    ----------
    table : hail.Table
        The table immediately after import and column standardization.
    file_path : str
        The source file path, for use in logging and error messages.
    validate_alleles : bool
        If True, validates that allele columns contain only ACGT characters.

    Returns
    -------
    hail.Table
        The fully processed and validated Hail Table.
    """
    if validate_alleles:
        table = _validate_alleles(table)
    table = _standardize_chromosome_column(table)
    table = table.filter(hl.is_defined(table.weight))
    filtered_row_count = table.count()
    if filtered_row_count == 0:
        raise ValueError(
            f"Input file '{file_path}' is empty or all variants were "
            f"filtered out due to missing weights or invalid alleles."
        )
    _check_duplicated_ids(table, file_path=file_path)
    logger.info(
        f"Successfully loaded {filtered_row_count} variants from {file_path}."
    )
    return table


def _read_prs_weights_noheader(
    file_path: str,
    column_map: dict,
    delimiter: str = '\t',
    keep_other_cols: bool = False,
    min_partitions: typing.Optional[int] = None,
    validate_alleles: bool = False
) -> hl.Table:
    """
    Reads a weight file without a header using a column map of indices.

    This internal function handles the specifics of importing a header-less
    file, mapping the provided integer indices to standardized column names,
    and then passing the resulting table to the main processing function.

    Parameters
    ----------
    file_path : str
        Path to the weight file.
    column_map : dict
        Dictionary mapping standard names to 1-based integer indices.
    delimiter : str
        Field delimiter.
    keep_other_cols : bool
        If True, all non-required columns are preserved.
    min_partitions : int, optional
        Hint for minimum number of partitions for the Hail Table.
    validate_alleles : bool
        If True, validates that allele columns contain only ACGT characters.

    Returns
    -------
    hail.Table
        The processed Hail Table.
    """
    logger.info(f"Importing file (no header): '{file_path}'")
    indices = list(column_map.values())
    if any(i < 1 for i in indices):
        raise ValueError(
            "Column indices in column_map must be 1-based and cannot be "
            "less than 1."
        )
    if len(indices) != len(set(indices)):
        raise ValueError("Duplicate column indices provided in column_map.")
    table = hl.import_table(
        file_path, delimiter=delimiter, no_header=True,
        min_partitions=min_partitions, comment='#'
    )
    standard_cols_exprs = {
        'chr': table[f"f{column_map['chr'] - 1}"],
        'pos': hl.int32(table[f"f{column_map['pos'] - 1}"]),
        'effect_allele': table[f"f{column_map['effect_allele'] - 1}"],
        'noneffect_allele': table[f"f{column_map['noneffect_allele'] - 1}"],
        'weight': hl.float64(table[f"f{column_map['weight'] - 1}"])
    }
    other_cols_exprs = {}
    if keep_other_cols:
        used_f_fields = {f'f{i - 1}' for i in indices}
        all_f_fields = table.row_value.dtype.fields
        other_fields = [f for f in all_f_fields if f not in used_f_fields]
        new_names = [f'non_req_col_{i+1}' for i, _ in enumerate(other_fields)]
        other_cols_exprs = {
            new: table[old] for new, old in zip(new_names, other_fields)
        }
    table = table.select(**standard_cols_exprs, **other_cols_exprs)
    return _process_prs_weights_table(table, file_path, validate_alleles)


def _read_prs_weights_header(
    file_path: str,
    column_map: dict,
    delimiter: str = '\t',
    keep_other_cols: bool = False,
    min_partitions: typing.Optional[int] = None,
    validate_alleles: bool = False
) -> hl.Table:
    """
    Reads a weight file with a header using a column map of names.

    This internal function handles the specifics of importing a file with a
    header, mapping the provided column names to standardized names, and
    then passing the resulting table to the main processing function.

    Parameters
    ----------
    file_path : str
        Path to the weight file.
    column_map : dict
        Dictionary mapping standard names to user-defined column names.
    delimiter : str
        Field delimiter.
    keep_other_cols : bool
        If True, all non-required columns are preserved.
    min_partitions : int, optional
        Hint for minimum number of partitions for the Hail Table.
    validate_alleles : bool
        If True, validates that allele columns contain only ACGT characters.

    Returns
    -------
    hail.Table
        The processed Hail Table.
    """
    logger.info(f"Importing file (with header): '{file_path}'")
    col_names = list(column_map.values())
    if len(col_names) != len(set(col_names)):
        raise ValueError("Duplicate column names provided in column_map.")
    types = {
        column_map['chr']: hl.tstr, column_map['pos']: hl.tint32,
        column_map['effect_allele']: hl.tstr,
        column_map['noneffect_allele']: hl.tstr,
        column_map['weight']: hl.tfloat64,
    }
    table = hl.import_table(
        file_path, delimiter=delimiter, no_header=False,
        min_partitions=min_partitions, types=types, comment='#', missing='NA'
    )
    missing = set(col_names) - set(table.row)
    if missing:
        raise ValueError(f"Required columns not in header: {missing}")
    standard_exprs = {
        'chr': table[column_map['chr']], 'pos': table[column_map['pos']],
        'effect_allele': table[column_map['effect_allele']],
        'noneffect_allele': table[column_map['noneffect_allele']],
        'weight': table[column_map['weight']]
    }
    other_exprs = {}
    if keep_other_cols:
        other_fields = [f for f in table.row if f not in col_names]
        other_exprs = {f: table[f] for f in other_fields}
    table = table.select(**standard_exprs, **other_exprs)
    return _process_prs_weights_table(table, file_path, validate_alleles)


def read_prs_weights(
    file_path: str,
    header: bool,
    column_map: dict,
    delimiter: str = '\t',
    keep_other_cols: bool = False,
    min_partitions: typing.Optional[int] = None,
    validate_alleles: bool = False
) -> hl.Table:
    """
    Reads a file containing variant effect weights for PRS calculation.

    This function requires a Hail-enabled environment to be active. It uses a
    flexible `column_map` dictionary to handle various input file formats.
    After standardizing the required columns, the function performs several
    validation checks, filtering out variants with missing weights, invalid
    alleles (if `validate_alleles=True`), or raising an error for duplicates.

    If a local file path is provided, it is automatically copied to a temporary
    directory in your GCS bucket for Hail to access.

    Parameters
    ----------
    file_path : str
        Path to the weight file (local or gs://).
    header : bool
        If True, `column_map` values should be strings (column names).
        If False, `column_map` values should be 1-based integers (column indices).
    column_map : dict
        A dictionary mapping standard names to user-defined names or indices.
        Must contain the keys: 'chr', 'pos', 'effect_allele',
        'noneffect_allele', and 'weight'.
        Example for header=True: {'chr': 'CHR', 'pos': 'BP', ...}
        Example for header=False: {'chr': 1, 'pos': 2, ...}
    delimiter : str, default '\t'
        Field delimiter.
    keep_other_cols: bool, default False
        If True, all columns not specified in `column_map` are preserved.
    min_partitions : int, optional
        A Hail-specific parameter to hint at the minimum number of partitions.
    validate_alleles : bool, default False
        If True, validates that allele columns contain only ACGT characters.

    Returns
    -------
    hail.Table
        A Hail Table with standardized columns.
    """
    gcs_path = _stage_local_file_to_gcs(file_path, sub_dir='temp_prs_data')
    required_keys = {
        'chr', 'pos', 'effect_allele', 'noneffect_allele', 'weight'
    }
    if not required_keys.issubset(column_map.keys()):
        missing = required_keys - set(column_map.keys())
        raise ValueError(f"column_map is missing required keys: {missing}")
    try:
        if hl.hadoop_stat(gcs_path)['size_bytes'] == 0:
            raise ValueError(f"Input file '{file_path}' is empty.")
    except hl.utils.java.FatalError as e:
        if 'Is a directory' not in str(e):
            raise e
    if header:
        if not all(isinstance(v, str) for v in column_map.values()):
            raise TypeError(
                "With header=True, column_map values must be strings."
            )
        return _read_prs_weights_header(
            file_path=gcs_path, column_map=column_map, delimiter=delimiter,
            keep_other_cols=keep_other_cols,
            min_partitions=min_partitions, validate_alleles=validate_alleles
        )
    else:
        if not all(isinstance(v, int) for v in column_map.values()):
            raise TypeError(
                "With header=False, column_map values must be integers."
            )
        return _read_prs_weights_noheader(
            file_path=gcs_path, column_map=column_map, delimiter=delimiter,
            keep_other_cols=keep_other_cols,
            min_partitions=min_partitions, validate_alleles=validate_alleles
        )


def read_prscs(file_path: str, **kwargs) -> hl.Table:
    """
    A simple wrapper to read PRS-CS output files.

    This function assumes a standard PRS-CS output format, which is a
    header-less, tab-separated file with the following columns:
    1. Chromosome
    2. Variant ID
    3. Base Position
    4. Effect Allele (A1)
    5. Non-Effect Allele (A2)
    6. Posterior Effect Size (weight)

    Note: The second column (Variant ID) is not loaded by default, as it is
    not required for the core functionality. To preserve this and any other
    columns, set `keep_other_cols=True` when calling this function.

    Parameters
    ----------
    file_path : str
        Path to the PRS-CS output file.
    **kwargs
        Other optional arguments to pass to `read_prs_weights`, such as
        `keep_other_cols` or `validate_alleles`.

    Returns
    -------
    hail.Table
        A processed Hail Table of the PRS-CS weights.
    """
    logger.info(f"Reading PRS-CS file: {file_path}")
    prscs_map = {
        'chr': 1, 'pos': 3, 'effect_allele': 4,
        'noneffect_allele': 5, 'weight': 6
    }
    return read_prs_weights(
        file_path=file_path, header=False, column_map=prscs_map, **kwargs
    )
