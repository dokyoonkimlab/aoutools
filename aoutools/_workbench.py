"""Environment helpers for the All of Us Researcher Workbench.

Two things about the Workbench are awkward enough, and change often enough,
that every notebook ends up re-deriving them. They live here so they exist in
exactly one place:

* Hail must be told which project pays for reads of the requester-pays bucket
  the VDS lives in, and the reference genome is no longer set through
  ``hl.init``.
* The environment variable that used to point at the VDS is not always set.

Neither function is required to use ``aoutools`` -- ``calculate_prs`` and
friends take a VDS you have already read. They exist to keep the boilerplate
short and correct.
"""

import logging
import os
import warnings

import hail as hl

logger = logging.getLogger(__name__)

# The path to the current All of Us short-read WGS VariantDataset.
#
# This is a *fallback*, and deliberately the last resort in `get_vds_path`. The
# Workbench used to export `WGS_VDS_PATH`; current images do not, so without a
# constant here every notebook hardcodes its own copy. As soon as the platform
# restores the variable, `get_vds_path` stops reading this and it becomes dead
# code -- which is the intent.
#
# It is pinned to a specific data release on purpose. Resolving "the newest
# version in the bucket" would remove the need to bump this by hand, but it
# could also hand back a genomics version that does not match the CDR the
# workspace is registered against -- a silently wrong analysis rather than an
# error. A stale constant fails loudly (the path does not exist); a wrong guess
# does not fail at all.
#
# Bump this when All of Us cuts a new data release.
DEFAULT_VDS_PATH = (
    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/vds/hail.vds"
)


def get_vds_path(path: str | None = None) -> str:
    """
    Resolves the path to the All of Us WGS VariantDataset.

    Checks three sources, in order, and returns the first one that is set:

    1. The `path` argument, if given.
    2. The `WGS_VDS_PATH` environment variable, exported by some Workbench
       images.
    3. `DEFAULT_VDS_PATH`, the current All of Us data release. Using this
       fallback emits a `UserWarning`, because it is pinned to a data version
       that this package cannot verify against your workspace.

    Parameters
    ----------
    path : str, optional
        An explicit path, which overrides both other sources. Passing the path
        you already have is always allowed, so a caller never has to branch on
        whether the environment happens to be configured.

    Returns
    -------
    str
        A GCS path to the VDS, suitable for `hl.vds.read_vds`.

    Examples
    --------
    >>> import hail as hl
    >>> from aoutools import get_vds_path, init_hail
    >>> init_hail()
    >>> vds = hl.vds.read_vds(get_vds_path())
    """
    if path:
        return path

    env_path = os.getenv("WGS_VDS_PATH")
    if env_path:
        return env_path

    warnings.warn(
        "The 'WGS_VDS_PATH' environment variable is not set, so aoutools is "
        f"falling back to its built-in default:\n  {DEFAULT_VDS_PATH}\n"
        "This is pinned to a specific All of Us data release. If your "
        "workspace uses a different release, pass the correct path explicitly "
        "-- an unnoticed mismatch between the genomic data and the CDR your "
        "phenotypes come from would produce a wrong analysis, not an error.",
        UserWarning,
        stacklevel=2,
    )
    return DEFAULT_VDS_PATH


def init_hail(reference: str = "GRCh38", **kwargs) -> None:
    """
    Initializes Hail for the All of Us Researcher Workbench.

    Wraps `hl.init` with the two settings the Workbench needs:

    * **Requester-pays billing.** The bucket holding the VDS charges the reader,
      so Hail must be given a project to bill. It is taken from the
      `GOOGLE_PROJECT` environment variable, which the Workbench sets. Without
      it, reading the VDS fails.
    * **The reference genome.** `hl.init(default_reference=...)` is deprecated;
      the reference is now set by calling `hl.default_reference` *with* an
      argument, after initialization. Doing it here keeps the deprecation
      warning out of your notebook.

    Calling this off the Workbench (where `GOOGLE_PROJECT` is unset) is fine:
    the requester-pays setting is simply omitted.

    This function is idempotent by default, so re-running the cell that calls it
    will not raise.

    Parameters
    ----------
    reference : str, default="GRCh38"
        The default reference genome. All of Us WGS data is GRCh38.
    **kwargs
        Passed through to `hl.init`. An explicit
        `gcs_requester_pays_configuration` or `idempotent` overrides the default
        set here.

    Examples
    --------
    >>> from aoutools import init_hail
    >>> init_hail()
    """
    kwargs.setdefault("idempotent", True)

    google_project = os.getenv("GOOGLE_PROJECT")
    if google_project:
        kwargs.setdefault("gcs_requester_pays_configuration", google_project)
    else:
        logger.warning(
            "The 'GOOGLE_PROJECT' environment variable is not set, so Hail is "
            "being initialized without a requester-pays billing project. "
            "Reading the All of Us VDS will fail. This is expected off the "
            "Researcher Workbench."
        )

    hl.init(**kwargs)
    # Must come after hl.init: passing default_reference to hl.init is
    # deprecated, and calling hl.default_reference() with no argument is the
    # getter. With an argument it is the setter.
    hl.default_reference(reference)
