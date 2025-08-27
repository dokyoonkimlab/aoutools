"""
PGS Catalog score file downloader using isolated pgscatalog.core CLI.

Supports local or GCS bucket output paths with temp staging for GCS.

This module internally manages a virtual environment to install
pgscatalog.core, avoiding dependency conflicts (tenacity version conflict
between pgscatalog.core and dsub).
"""

import concurrent.futures
import importlib.util
import logging
import os
import pathlib
import subprocess
import sys
import tempfile
import threading
import venv
from pathlib import Path
from typing import Iterable, Optional, Union

from google.cloud import storage
from google.cloud.exceptions import GoogleCloudError

logger = logging.getLogger(__name__)

_cached_cli_path_lock = threading.Lock()
_cached_cli_path: Optional[Path] = None

DEFAULT_ENV_DIR = Path.home() / ".aoutools" / "pgscatalog_env"
PGS_ENV_DIR = Path(os.environ.get("AOUTOOLS_PGS_ENV_DIR", DEFAULT_ENV_DIR))


def _run(cmd: list[str], **kwargs) -> None:
    """Run a shell command and raise an error if it fails, logging output."""
    try:
        completed = subprocess.run(
            cmd, check=True, capture_output=True, text=True, **kwargs
        )
        logger.debug("Command output: %s", completed.stdout)
    except subprocess.CalledProcessError as e:
        logger.error("Command '%s' failed with exit code %d", cmd, e.returncode)
        logger.error("stdout:\n%s", e.stdout)
        logger.error("stderr:\n%s", e.stderr)
        raise


def _create_env() -> Path:
    """Create an empty virtual environment and return its bin directory."""
    PGS_ENV_DIR.parent.mkdir(parents=True, exist_ok=True)
    venv.create(PGS_ENV_DIR, with_pip=True)
    return _get_bin_dir()


def _get_bin_dir() -> Path:
    """Return the path to the bin/Scripts directory of the isolated venv."""
    return PGS_ENV_DIR / ("Scripts" if sys.platform == "win32" else "bin")


def _get_pgscatalog_version(bin_dir: Path) -> Optional[str]:
    """
    Get the installed pgscatalog.core version from the isolated venv.

    Returns None if not installed or if query fails.
    """
    try:
        output = subprocess.check_output(
            [bin_dir / "python", "-m", "pip", "show", "pgscatalog.core"],
            text=True,
        )
        for line in output.splitlines():
            if line.startswith("Version:"):
                return line.split(":", 1)[1].strip()
    except subprocess.CalledProcessError:
        return None
    return None


def _ensure_pgscatalog_download(min_version: str = "1.0.1") -> Path:
    """
    Ensure pgscatalog.core is installed in an isolated venv.

    If the environment does not exist, it will be created.
    If the installed version is older than min_version, it will be upgraded.
    If the environment is corrupted, it will be recreated.

    Returns
    -------
    Path
        Path to the pgscatalog-download CLI executable inside the venv.
    """
    global _cached_cli_path
    with _cached_cli_path_lock:
        if _cached_cli_path is not None:
            return _cached_cli_path

        bin_dir = _get_bin_dir()
        cli_path = bin_dir / "pgscatalog-download"

        if not bin_dir.exists() or not (bin_dir / "pip").exists():
            logger.info(
                "Creating isolated environment for pgscatalog.core >= %s",
                min_version
            )
            bin_dir = _create_env()

        installed_version = _get_pgscatalog_version(bin_dir)

        if installed_version is None:
            logger.info("Installing pgscatalog.core >= %s", min_version)
            _run([
                bin_dir / "pip",
                "install",
                "--no-user",
                f"pgscatalog.core>={min_version}"
            ])
        elif installed_version < min_version:
            logger.info(
                "Upgrading pgscatalog.core to >= %s (current: %s)",
                min_version,
                installed_version,
            )
            _run([
                bin_dir / "pip",
                "install",
                "--no-user",
                "--upgrade",
                f"pgscatalog.core>={min_version}"
            ])

        _cached_cli_path = cli_path
        return cli_path


def _run_pgscatalog_download(
    outdir: Union[str, Path],
    *args: str,
    min_version: str = "1.0.1",
) -> None:
    """
    Run a pgscatalog.core CLI command in its isolated environment.

    Parameters
    ----------
    outdir : str or Path
        Mandatory output directory path for pgscatalog-download `-o` option.
    *args : str
        Other CLI arguments to pass to pgscatalog.
    min_version : str, optional
        Minimum version of pgscatalog.core to ensure.
    """
    cli_path = _ensure_pgscatalog_download(min_version=min_version)
    outdir_path = str(outdir)
    cmd = [str(cli_path), "-o", outdir_path, *args]
    _run(cmd)


def _normalize_arg(arg: Union[Iterable[str], str, None]) -> list[str]:
    if arg is None:
        return []
    if isinstance(arg, str):
        return [arg]
    return list(arg)


def download_pgs(
    *,
    outdir: Union[str, pathlib.Path],
    pgs: Union[Iterable[str], str, None] = None,
    efo: Union[Iterable[str], str, None] = None,
    pgp: Union[Iterable[str], str, None] = None,
    build: Optional[str] = "GRCh38",
    efo_include_children: bool = True,
    overwrite_existing_file: bool = False,
    user_agent: Optional[str] = None,
    verbose: bool = False,
) -> None:
    """
    Download PGS Catalog scoring files to a local directory or GCS bucket.

    This function detects if the output path is local or GCS (gs://).
    If GCS, it downloads files to a temp directory then uploads to GCS.

    Parameters
    ----------
    outdir : str or pathlib.Path
        Local directory or GCS bucket path (e.g., 'gs://my-bucket/path').
    pgs : str or iterable of str, optional
        PGS Catalog ID(s) (e.g., "PGS000194").
    efo : str or iterable of str, optional
        EFO term(s) (e.g., "EFO_0004611").
    pgp : str or iterable of str, optional
        PGP publication ID(s).
    build : str, optional
        Genome build ("GRCh37" or "GRCh38"), default "GRCh38".
    efo_include_children : bool, default True
        Whether to include descendant EFO terms.
    overwrite_existing_file : bool, default False
        Overwrite existing files if newer versions exist.
    user_agent : str, optional
        Custom user agent string.
    verbose : bool, default False
        Enable verbose logging.

    Returns
    -------
    None
        The PGS Catalog score file(s) saved to the specified output path.

    Raises
    ------
    FileNotFoundError
        If local output directory does not exist.
    ValueError
        If none of pgs, efo, or pgp are provided.
    Exception
        On download or upload failure.
    """
    pgs_args = _normalize_arg(pgs)
    efo_args = _normalize_arg(efo)
    pgp_args = _normalize_arg(pgp)

    if not (pgs_args or efo_args or pgp_args):
        raise ValueError(
            "At least one of 'pgs', 'efo', or 'pgp' must be provided."
        )

    cli_args = []
    if pgs_args:
        cli_args.extend(["-i", *pgs_args])
    if efo_args:
        cli_args.extend(["-t", *efo_args])
    if pgp_args:
        cli_args.extend(["-p", *pgp_args])
    if build:
        cli_args.extend(["-b", build])
    if efo_include_children:
        cli_args.append("-e")
    if overwrite_existing_file:
        cli_args.append("-w")
    if user_agent:
        cli_args.extend(["-c", user_agent])

    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("pgscatalog.corelib").setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    outdir_str = str(outdir)

    if outdir_str.startswith("gs://"):
        with tempfile.TemporaryDirectory() as temp_dir:
            logger.info(
                "Downloading scoring files to temporary directory: %s",
                temp_dir
            )
            _run_pgscatalog_download(temp_dir, *cli_args)

            storage_client = storage.Client()
            bucket_name, *path_parts = outdir_str[5:].split("/", 1)
            prefix = path_parts[0] if path_parts else ""
            bucket = storage_client.bucket(bucket_name)

            logger.info("Uploading files to gs://%s/%s", bucket_name, prefix)

            files_to_upload = list(pathlib.Path(temp_dir).iterdir())
            with concurrent.futures.ThreadPoolExecutor() as executor:
                futures = []
                for local_path in files_to_upload:
                    blob_name = os.path.join(prefix, local_path.name)
                    blob = bucket.blob(blob_name)
                    future = executor.submit(
                        blob.upload_from_filename, str(local_path)
                    )
                    futures.append(future)

                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()  # Raises exceptions from the thread
                    except GoogleCloudError as e:
                        logger.error("Failed to upload file to GCS: %s", e)
                        executor.shutdown(wait=False, cancel_futures=True)
                        raise

            logger.info(
                "Upload of %d files to GCS completed.", len(files_to_upload)
            )

    else:
        # Local path: verify dir exists, then download directly
        local_path = pathlib.Path(outdir_str).expanduser()
        if not local_path.exists() or not local_path.is_dir():
            error_msg = (
                f"Local output directory '{local_path}' does not exist or "
                "is not a directory."
            )
            raise FileNotFoundError(error_msg)

        logger.info(
            "Downloading scoring files directly to local directory: %s",
            local_path
        )
        _run_pgscatalog_download(local_path, *cli_args)
        logger.info("Download to local directory complete.")
