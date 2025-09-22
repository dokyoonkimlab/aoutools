"""
High-level workflow function for end-to-end PGS calculation.
"""

import logging
import tempfile
from pathlib import Path
from typing import Optional, Union, Iterable

import hail as hl

from ._downloader import download_pgs
from ._reader import read_prs_weights
from ._calculator_batch import calculate_prs_batch
from ._config import PRSConfig

logger = logging.getLogger(__name__)


def calculate_pgs(
    *,
    vds: hl.vds.VariantDataset,
    output_path: str,
    pgs: Union[Iterable[str], str],
    build: Optional[str] = "GRCh38",
    config: PRSConfig = PRSConfig(),
    user_agent: Optional[str] = None,
    verbose: bool = False,
) -> Optional[str]:
    """
    Downloads specified PGS Catalog scoring files and calculates PRS.

    This function automates a controlled workflow:

    1. **Download**: Fetches scoring files from the PGS Catalog for a specific
        list of PGS IDs.
    2.  **Read**: Parses the downloaded scoring files into Hail Tables.
    3. **Calculate**: Computes the Polygenic Risk Score(s) for each downloaded
        file and exports a single CSV file.

    Notes
    -----
    This function does not accept EFO traits or PGP publication IDs to prevent
    the unexpected download of a large number of scoring files, which may
    overwhelm the storage or computational resources of a typical workspace.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        A Hail VariantDataset containing the genotype data to be scored.
    output_path : str
        A GCS path (e.g., 'gs://bucket/results.csv') for the output file.
    pgs : str or iterable of str
        One or more PGS Catalog ID(s) (e.g., "PGS000771") to download.
        This argument is required.
    build : str, optional
        The genome build for harmonized scores ("GRCh37" or "GRCh38").
        Defaults to "GRCh38".
    config : PRSConfig, optional
        A configuration object for calculation parameters.
    user_agent : str, optional
        A custom user agent string for PGS Catalog API requests.
    verbose : bool, default False
        Enable verbose logging for the download process.

    Returns
    -------
    str or None
        The output path if results are successfully written; otherwise, None.

    Raises
    ------
    ValueError
        If the `output_path` is not a valid GCS path, or if a downloaded
        scoring file is empty, malformed, or contains duplicate variants.
    TypeError
        If `config.samples_to_keep` is an unsupported type.
    Exception
        If the download process fails due to network issues, invalid PGS IDs,
        or other errors from the underlying `pgscatalog-download` tool.
    """
    # Define the column mapping for standard PGS Catalog scoring files
    pgs_column_map = {
        'chr': 'hm_chr',
        'pos': 'hm_pos',
        'effect_allele': 'effect_allele',
        'noneffect_allele': 'other_allele',
        'weight': 'effect_weight',
    }

    with tempfile.TemporaryDirectory() as temp_dir:
        logger.info(
            "Created temporary directory for PGS downloads: %s", temp_dir
        )

        # Step 1: Download scoring files for specified PGS IDs
        download_pgs(
            outdir=temp_dir,
            pgs=pgs,
            build=build,
            user_agent=user_agent,
            verbose=verbose,
            overwrite_existing_file=True,
        )

        # Step 2: Read downloaded files and prepare for batch calculation
        downloaded_files = list(Path(temp_dir).rglob("*.txt.gz"))
        if not downloaded_files:
            logger.warning(
                "No scoring files were downloaded for the specified PGS IDs."
            )
            return None

        weights_tables_map = {}
        for file_path in downloaded_files:
            score_name = file_path.name.split('.')[0]
            logger.info(
                "Reading scoring file for '%s' from %s", score_name, file_path
            )
            try:
                weights_table = read_prs_weights(
                    file_path=str(file_path),
                    header=True,
                    column_map=pgs_column_map,
                    delimiter='\t',
                    comment='#',
                    validate_alleles=True,
                    missing='',
                    force=True,
                )
                weights_tables_map[score_name] = weights_table
            except (ValueError, TypeError) as e:
                logger.warning(
                    "Could not read scoring file for '%s'. Skipping. Error: %s. "
                    "Please check the format of the downloaded file; some "
                    "scoring files may contain variants without positional "
                    "information (e.g., HLA types) or with incomplete effect "
                    "allele information, which are not supported.",
                    score_name, e
                )

        if not weights_tables_map:
            logger.warning("No valid scoring files could be read. Aborting.")
            return None

        # Step 3: Calculate PRS for all scores in batch mode
        logger.info(
            "Calculating PRS for %d scores in batch mode.",
            len(weights_tables_map)
        )
        result_path = calculate_prs_batch(
            weights_tables_map=weights_tables_map,
            vds=vds,
            output_path=output_path,
            config=config,
        )

        return result_path
