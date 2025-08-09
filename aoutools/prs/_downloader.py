"""Downloader for PGS Catalog scoring files."""

from typing import Iterable, Union
import concurrent
import functools
import logging
import pathlib
import time
from concurrent.futures import ThreadPoolExecutor

import requests
from tqdm import tqdm

from pgscatalog.core.lib import ScoringFiles, GenomeBuild, Config

logger = logging.getLogger(__name__)


def _normalize_arg(arg: Union[Iterable[str], str, None]) -> list[str]:
    if arg is None:
        return []
    if isinstance(arg, str):
        return [arg]
    if isinstance(arg, Iterable):
        return list(arg)
    raise TypeError(f"Unsupported type for argument: {type(arg)}")


def _rate_limited(max_calls_per_minute):
    """Decorator to limit function calls to max_calls_per_minute."""
    min_interval = 60.0 / max_calls_per_minute

    def decorator(func):
        last_called = [0.0]

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            elapsed = time.time() - last_called[0]
            wait = min_interval - elapsed
            if wait > 0:
                logger.info("Rate limit reached, sleeping %.2f seconds", wait)
                time.sleep(wait)
            result = func(*args, **kwargs)
            last_called[0] = time.time()
            return result

        return wrapper

    return decorator



def _retry_on_429(max_retries=3, initial_delay=1):
    """Decorator to retry a function if HTTP 429 is raised."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            delay = initial_delay
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except requests.HTTPError as e:
                    if e.response.status_code == 429:
                        logger.warning(
                            "Rate limit hit (HTTP 429), retrying after %s seconds",
                            delay,
                        )
                        time.sleep(delay)
                        delay *= 2
                    else:
                        raise
            # Last try (raises if it fails)
            return func(*args, **kwargs)

        return wrapper

    return decorator


@_rate_limited(100)  # max 100 calls per minute to the API
@_retry_on_429(max_retries=5, initial_delay=1)
def fetch_scoring_files(*args, **kwargs):
    """
    Wrapper around ScoringFiles constructor to apply rate limiting
    and retry on 429.
    """
    return ScoringFiles(*args, **kwargs)


def download_pgs(
    outdir: str,
    pgs: Union[Iterable[str], str, None] = None,
    efo: Union[Iterable[str], str, None] = None,
    pgp: Union[Iterable[str], str, None] = None,
    build: str | None = "GRCh38",
    efo_include_children: bool = True,
    overwrite_existing_file: bool = False,
    user_agent: str | None = None,
    verbose: bool = False,
) -> None:
    """
    Download PGS Catalog scoring files asynchronously.

    Parameters
    ----------
    outdir : str
        A directory path where downloaded scoring files will be saved.
    pgs : str or iterable of str, optional
        PGS Catalog ID(s) (e.g., "PGS000194"). Can be a single string, list,
        tuple, or set of strings.
    efo : str or iterable of str, optional
        Traits described by EFO term(s) (e.g., "EFO_0004611"). Can be a single
        string or iterable of strings.
    pgp : str or iterable of str, optional
        PGP publication ID(s) (e.g., "PGP000007"). Can be a single string or
        iterable of strings.
    build : str, optional
        Genome build for harmonized scores: "GRCh37" or "GRCh38". Default is
        "GRCh38". Choosing "GRCh37" triggers a warning as All of Us genetic
        data is based on GRCh38.
    efo_include_children : bool, default True
        Whether to include scoring files tagged with descendant EFO terms.
    overwrite_existing_file : bool, default False
        Whether to overwrite existing files if newer versions are available.
    user_agent : str, optional
        A custom user agent string for PGS Catalog API requests.
    verbose : bool, default False
        Enable verbose logging output.

    Returns
    -------
    None
        This function does not return anything.

    Raises
    ------
    FileNotFoundError
        If the output directory does not exist.
    ValueError
        If none of `pgs`, `efo`, or `pgp` parameters are provided.
    Exception
        If any download task raises an exception.
    """
    pgs = _normalize_arg(pgs)
    efo = _normalize_arg(efo)
    pgp = _normalize_arg(pgp)

    if not pgs and not efo and not pgp:
        raise ValueError(
            "At least one of pgs, efo, or pgp must be provided"
        )

    if verbose:
        logging.getLogger("pgscatalog.corelib").setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    outdir_path = pathlib.Path(outdir).expanduser()
    if not outdir_path.exists():
        raise FileNotFoundError(
            "Output directory '%s' does not exist", outdir
        )

    if user_agent is not None:
        logger.info("Setting user agent to %s", user_agent)
        Config.API_HEADER = {"user-agent": user_agent}

    if build == "GRCh37":
        logger.warning(
            "Warning: All of Us Genetic data is based on GRCh38, "
            "but GRCh37 was specified."
        )

    build_enum = GenomeBuild.from_string(build) if build else None
    if build_enum is None and build is not None:
        logger.warning(
            "Invalid genome build '%s', proceeding without harmonized build",
            build,
        )
    else:
        logger.info(
            "Downloading scoring files harmonized to build: %s", build_enum
        )

    # Use rate-limited + retry wrapper here
    sfs = fetch_scoring_files(
        [*pgs, *pgp, *efo],
        target_build=build_enum,
        include_children=efo_include_children,
    )

    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for scorefile in sfs:
            logger.info("Submitting %r download", scorefile)
            futures.append(
                executor.submit(
                    scorefile.download,
                    overwrite=overwrite_existing_file,
                    directory=outdir_path,
                )
            )

        for future in tqdm(
            concurrent.futures.as_completed(futures), total=len(futures)
        ):
            future.result()
            logger.info("Download complete")

    logger.info("All downloads finished")
