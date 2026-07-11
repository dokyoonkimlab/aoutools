"""Utility functions for Hail data processing"""

import logging
import os
from contextlib import contextmanager
from time import perf_counter

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
        raise OSError(
            "The 'WORKSPACE_BUCKET' environment variable is not set. "
            "This is required to stage local files to GCS."
        )

    gcs_path = os.path.join(
        workspace_bucket, 'data', sub_dir, os.path.basename(file_path)
    )

    logger.info(
        "Local file detected. Staging '%s' to '%s'...",
        file_path,
        gcs_path,
    )
    hfs.copy(f'file://{os.path.abspath(file_path)}', gcs_path)

    return gcs_path


def _standardize_chromosome_column(table: hl.Table) -> hl.Table:
    """
    Ensures that every value in the 'chr' column has a 'chr' prefix.

    Prepends the 'chr' prefix to any value that lacks it (e.g., '1' becomes
    'chr1'), while leaving already-prefixed values (e.g., 'chr1') unchanged.
    The check is applied per row rather than inferred from a single sampled
    value, so mixed-format inputs (some rows prefixed, some not) are handled
    correctly. This standardization is crucial for matching against reference
    datasets like the All of Us VDS.

    Parameters
    ----------
    table : hail.Table
        A Hail Table to process; must contain a 'chr' column.

    Returns
    -------
    hail.Table
        A Hail Table with a standardized 'chr' column.
    """
    # Cast to string so integer contigs (e.g. 1 rather than '1') are tolerated.
    chr_str = hl.str(table.chr)
    # Per-row and idempotent: already-prefixed values pass through untouched,
    # so this is safe on mixed-format columns and if applied more than once.
    return table.annotate(
        chr=hl.if_else(chr_str.startswith('chr'), chr_str, 'chr' + chr_str)
    )
