"""Utility functions for Hail data processing"""

import os
import logging
from time import perf_counter
from contextlib import contextmanager
import hail as hl
import hailtop.fs as hfs

# Configure a logger for module-level use.
logger = logging.getLogger(__name__)


@contextmanager
def _log_timing(description: str, enabled: bool = True):
    """
    Context manager to log the duration of a code block.

    Parameters
    ----------
    description : str
        A description of the code block being timed.
    enabled : bool, default=True
        If False, disables timing and logging.
    """
    if not enabled:
        yield
        return

    start_time = perf_counter()
    logger.info("%s...", description)
    yield
    duration = perf_counter() - start_time
    logger.info("%s finished in %.2f seconds.", description, duration)


def _stage_local_file_to_gcs(
    file_path: str,
    sub_dir: str
) -> str:
    """
    Checks if file path is local; if so, stages it to GCS.

    Useful for distributed platforms like All of Us Researcher Workbench, where
    Hail's Spark cluster cannot access local notebook environment directly. This
    function copies local files to a subdirectory within workspace's GCS bucket,
    which is accessible by the cluster.

    Uses WORKSPACE_BUCKET environment variable provided by the platform.

    Parameters
    ----------
    file_path : str
        A path to the file.
    sub_dir : str
        A subdirectory within `$WORKSPACE_BUCKET/data/` to copy the file to.

    Returns
    -------
    str
        A GCS path to the file accessible by Hail.

    Raises
    ------
    FileNotFoundError
        If a local `file_path` does not exist.
    EnvironmentError
        If the `WORKSPACE_BUCKET` environment variable is not set.
    """
    if file_path.startswith('gs://'):
        return file_path
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Local file does not exist: {file_path}")

    workspace_bucket = os.getenv('WORKSPACE_BUCKET')
    if not workspace_bucket:
        raise EnvironmentError(
            "The 'WORKSPACE_BUCKET' environment variable is not set. "
            "This is required to stage local files to GCS."
        )

    gcs_path = os.path.join(
        workspace_bucket, 'data', sub_dir, os.path.basename(file_path)
    )

    if hfs.exists(gcs_path):
        logger.info("File already exists at %s, skipping copy.", gcs_path)
    else:
        logger.info(
            "Local file detected. Staging '%s' to '%s'...",
            file_path,
            gcs_path,
        )
        hfs.copy(f'file://{os.path.abspath(file_path)}', gcs_path)

    return gcs_path


def _standardize_chromosome_column(table: hl.Table) -> hl.Table:
    """
    Ensures that the 'chr' column has a 'chr' prefix.

    Inspects a sample value from the 'chr' column. If the prefix is missing
    (e.g., '1' instead of 'chr1'), annotates the entire column to add it.
    This standardization is crucial for matching against reference datasets
    like the All of Us VDS.

    Parameters
    ----------
    table : hail.Table
        A Hail Table to process; must contain a 'chr' column.

    Returns
    -------
    hail.Table
        A Hail Table with a standardized 'chr' column.
    """
    if table.count() == 0:
        return table

    sample_chr = table.select('chr').take(1)[0].chr
    if not str(sample_chr).startswith('chr'):
        logger.info("Adding 'chr' prefix to chromosome column.")
        table = table.annotate(chr=hl.str('chr') + table.chr)

    return table
