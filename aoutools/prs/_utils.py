"""Utility functions for Hail data processing"""

import os
import typing
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
    A context manager to log the duration of a code block.

    Parameters
    ----------
    description : str
        A description of the code block being timed.
    enabled : bool
        If False, the timer is disabled.
    """
    if not enabled:
        yield
        return

    start_time = perf_counter()
    logger.info("%s...", description)
    yield
    duration = perf_counter() - start_time
    logger.info("%s finished in %.2f seconds.", description, duration)


def _stage_local_file_to_gcs(file_path: str, sub_dir: str) -> str:
    """
    Checks if a file path is local and, if so, stages it to GCS.

    For distributed platforms like the All of Us Researcher Workbench, Hail's
    Spark cluster cannot directly access the local notebook environment. This
    function ensures that local files are copied to a subdirectory within
    the workspace's GCS bucket, which is accessible by the cluster.

    It uses the WORKSPACE_BUCKET environment variable provided by the platform.

    Parameters
    ----------
    file_path : str
        The path to the file.
    sub_dir : str
        The subdirectory within `$WORKSPACE_BUCKET/data/` to copy the file to.

    Returns
    -------
    str
        A GCS path to the file that Hail can access.

    Raises
    ------
    FileNotFoundError
        If a local `file_path` is provided and the file does not exist.
    EnvironmentError
        If a local file is provided and the `WORKSPACE_BUCKET` environment
        variable is not set.
    """
    if file_path.startswith('gs://'):
        return file_path
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Local file does not exist: {file_path}")

    workspace_bucket = os.getenv('WORKSPACE_BUCKET')
    # Fail fast if the required environment variable is not set.
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
    Ensures the 'chr' column has a 'chr' prefix.

    This function inspects a sample of the 'chr' column. If the 'chr' prefix
    is missing (e.g., '1' instead of 'chr1'), it annotates the entire
    column to add the prefix. This is crucial for matching against reference
    datasets like the All of Us VDS.

    Parameters
    ----------
    table : hail.Table
        The table to process, must contain a 'chr' column.

    Returns
    -------
    hail.Table
        A table with a standardized 'chr' column.
    """
    if table.count() == 0:
        return table

    sample_chr = table.select('chr').take(1)[0].chr
    if not str(sample_chr).startswith('chr'):
        logger.info("Adding 'chr' prefix to chromosome column.")
        table = table.annotate(chr=hl.str('chr') + table.chr)

    return table


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
