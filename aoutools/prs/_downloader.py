"""
PGS Catalog score file downloader using isolated pgscatalog.core CLI.

Supports local or GCS bucket output paths with temp staging for GCS.

This module internally manages a virtual environment to install
pgscatalog.core, avoiding dependency conflicts (tenacity version conflict
between pgscatalog.core and dsub).
"""

import concurrent.futures
import logging
import os
import pathlib
import subprocess
import sys
import tempfile
import threading
import venv
from collections.abc import Iterable
from pathlib import Path

from google.cloud import storage
from google.cloud.exceptions import GoogleCloudError
from packaging.version import Version

logger = logging.getLogger(__name__)

_cached_cli_path_lock = threading.Lock()
_cached_cli_path: Path | None = None

DEFAULT_ENV_DIR = Path.home() / ".aoutools" / "pgscatalog_env"
PGS_ENV_DIR = Path(os.environ.get("AOUTOOLS_PGS_ENV_DIR", DEFAULT_ENV_DIR))


def _subprocess_env() -> dict[str, str]:
    """
    A copy of the environment scrubbed of settings that defeat the venv.

    The isolated venv gives pgscatalog.core its own site-packages, but that
    only isolates *imports*. Two settings the Workbench image can inherit still
    let the base environment's pinned `tenacity` (8.2.3, held by `dsub`)
    override the `>=9.0.0` that `pgscatalog.core` needs:

    - `PIP_CONSTRAINT` makes pip apply the base image's constraints file inside
      the venv, so `pip install pgscatalog.core` fails to resolve tenacity.
    - `PYTHONPATH` prepends the base site-packages onto `sys.path`, so the CLI
      imports the base tenacity at runtime even after a clean install.

    Stripping those -- plus `PYTHONHOME`/`PIP_USER`, and forcing
    `PYTHONNOUSERSITE` -- makes the child processes behave like pipx's: they see
    only the venv. Index configuration (`PIP_INDEX_URL` and friends) is left
    alone, since a Workbench mirror is load-bearing for the download itself.
    """
    env = os.environ.copy()
    for var in ("PIP_CONSTRAINT", "PIP_USER", "PYTHONPATH", "PYTHONHOME"):
        env.pop(var, None)
    env["PYTHONNOUSERSITE"] = "1"
    return env


def _run(cmd: list[str], **kwargs) -> None:
    """
    Run a shell command, raising with the subprocess's own error text.

    The failure text is put into the **exception message**, not only into the
    log. `aoutools` installs a `NullHandler`, and a notebook configures no
    logging by default, so anything logged here is discarded -- which used to
    leave a caller with a bare `CalledProcessError` and no way at all to see
    why the download failed.

    Runs with a scrubbed environment (see `_subprocess_env`) unless the caller
    overrides `env`, so the isolated venv is not undone by an inherited
    `PIP_CONSTRAINT` or `PYTHONPATH`.
    """
    kwargs.setdefault("env", _subprocess_env())
    try:
        completed = subprocess.run(
            cmd, check=True, capture_output=True, text=True, **kwargs
        )
        logger.debug("Command output: %s", completed.stdout)
    except subprocess.CalledProcessError as e:
        # No logger.error here: the NullHandler default swallows it, and the
        # full failure text (exit code, stdout, stderr) goes into the raised
        # RuntimeError below, which is where a caller actually sees it.
        details = "\n".join(
            part
            for part in (
                f"Command failed with exit code {e.returncode}:",
                f"  {' '.join(str(c) for c in cmd)}",
                f"\nstderr:\n{e.stderr.strip()}" if e.stderr else "",
                f"\nstdout:\n{e.stdout.strip()}" if e.stdout else "",
            )
            if part
        )
        raise RuntimeError(details) from e


def _create_env() -> Path:
    """Create an empty virtual environment and return its bin directory."""
    PGS_ENV_DIR.parent.mkdir(parents=True, exist_ok=True)
    venv.create(PGS_ENV_DIR, with_pip=True)
    return _get_bin_dir()


def _get_bin_dir() -> Path:
    """Return the path to the bin/Scripts directory of the isolated venv."""
    return PGS_ENV_DIR / ("Scripts" if sys.platform == "win32" else "bin")


def _get_pgscatalog_version(bin_dir: Path) -> str | None:
    """
    Get the installed pgscatalog.core version from the isolated venv.

    Returns None if not installed or if query fails.
    """
    try:
        output = subprocess.check_output(
            [bin_dir / "python", "-m", "pip", "show", "pgscatalog.core"],
            text=True,
            env=_subprocess_env(),
        )
        for line in output.splitlines():
            if line.startswith("Version:"):
                return line.split(":", 1)[1].strip()
    except subprocess.CalledProcessError:
        return None
    return None


def _pgscatalog_constraint(min_version: str) -> str:
    """Return the version constraint string for pgscatalog.core."""
    if sys.version_info < (3, 11):
        return f"pgscatalog.core>={min_version},<1.0.2"
    return f"pgscatalog.core>={min_version}"


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
                min_version,
            )
            bin_dir = _create_env()

        installed_version = _get_pgscatalog_version(bin_dir)

        if installed_version is None:
            logger.info("Installing pgscatalog.core >= %s", min_version)
            _run(
                [
                    bin_dir / "pip",
                    "install",
                    "--no-user",
                    _pgscatalog_constraint(min_version),
                ]
            )
        elif Version(installed_version) < Version(min_version):
            logger.info(
                "Upgrading pgscatalog.core to >= %s (current: %s)",
                min_version,
                installed_version,
            )
            _run(
                [
                    bin_dir / "pip",
                    "install",
                    "--no-user",
                    "--upgrade",
                    _pgscatalog_constraint(min_version),
                ]
            )

        # Verify the CLI is actually there before caching it. Neither the venv
        # creation nor the pip install is checked otherwise, so a partial or
        # failed install would be cached as success and surface later as an
        # inscrutable error from a command that does not exist -- or worse, from
        # a stale one left behind by an earlier attempt.
        if not cli_path.exists():
            raise RuntimeError(
                f"pgscatalog.core was installed into {PGS_ENV_DIR}, but the "
                f"'pgscatalog-download' command is not there. The environment "
                f"is likely half-built from an interrupted install.\n\n"
                f"Delete it and let aoutools rebuild it:\n"
                f"  import shutil; shutil.rmtree('{PGS_ENV_DIR}')"
            )

        _cached_cli_path = cli_path
        return cli_path


def _run_pgscatalog_download(
    outdir: str | Path,
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

    try:
        _run(cmd)
    except RuntimeError as e:
        # Without `-w`, pgscatalog-download raises FileExistsError on the first
        # file that is already present -- and because the downloads run
        # concurrently, that kills the WHOLE batch, including the scores that
        # would have downloaded fine. Re-running a notebook cell is the normal
        # way to hit this, and the raw traceback (deep inside tenacity and a
        # thread pool) says nothing about how to fix it.
        if "already exists" not in str(e):
            raise
        raise FileExistsError(
            f"A scoring file already exists in {outdir_path}, and "
            "`pgscatalog-download` refuses to overwrite it -- which aborts "
            "the entire batch, not just that one score.\n\n"
            "Either pass `overwrite_existing_file=True` to re-download it, or "
            "delete the directory and start clean.\n\n"
            f"Original error:\n{e}"
        ) from e


def _normalize_arg(arg: Iterable[str] | str | None) -> list[str]:
    if arg is None:
        return []
    if isinstance(arg, str):
        return [arg]
    return list(arg)


def download_pgs(
    *,
    outdir: str | pathlib.Path,
    pgs: Iterable[str] | str | None = None,
    efo: Iterable[str] | str | None = None,
    pgp: Iterable[str] | str | None = None,
    build: str | None = "GRCh38",
    efo_include_children: bool = True,
    overwrite_existing_file: bool = True,
    user_agent: str | None = None,
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
    overwrite_existing_file : bool, default True
        Re-download a scoring file that is already present in `outdir`.

        With the default (True), re-running a cell that already downloaded its
        files simply fetches them again and replaces them. This is safe: the
        file for a given PGS Catalog ID never changes, so the new copy is
        identical to the old one. Pass False to make an already-present file an
        **error** instead -- but note that because the downloads run
        concurrently, that error aborts the whole batch, not just the score
        whose file was present.
    user_agent : str, optional
        Custom user agent string.

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

    Notes
    -----
    The first call in a session provisions an isolated Python virtual
    environment and ``pip install``\\ s ``pgscatalog.core`` into it, to avoid a
    ``tenacity`` version conflict between ``pgscatalog.core`` and ``dsub``.
    Consequently the first call:

    - requires **network access** and a **writable home directory**, and
    - incurs a one-time setup delay while the environment is built.

    Subsequent calls reuse the cached environment and skip this step. The
    environment lives at ``~/.aoutools/pgscatalog_env`` by default; set the
    ``AOUTOOLS_PGS_ENV_DIR`` environment variable (before the first call) to
    relocate it. To force a rebuild, delete that directory. See the
    "Using the ``download_pgs`` function" how-to guide for details.
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

    outdir_str = str(outdir)

    if outdir_str.startswith("gs://"):
        with tempfile.TemporaryDirectory() as temp_dir:
            logger.debug(
                "Downloading scoring files to temporary directory: %s", temp_dir
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
                    except GoogleCloudError:
                        # Let the exception propagate as-is; it already names
                        # the failed upload. A logger.error here would only
                        # double-report under a configured handler and vanish
                        # under the NullHandler default.
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
            local_path,
        )
        _run_pgscatalog_download(local_path, *cli_args)
        logger.info("Download to local directory complete.")
