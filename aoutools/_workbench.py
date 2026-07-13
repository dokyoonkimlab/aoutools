"""Environment helpers for the All of Us Researcher Workbench.

A few things about the Workbench are awkward enough, and change often enough,
that every notebook ends up re-deriving them. They live here so they exist in
exactly one place:

* Hail must be told which project pays for reads of the requester-pays bucket
  the VDS lives in, and the reference genome is no longer set through
  ``hl.init``.
* The environment variables that used to point at the VDS and at the workspace
  bucket are **not set at all** on the current (Verily) platform. The older
  Workbench guaranteed them; code that assumed they exist now raises `KeyError`
  or `OSError` on a fresh workspace.

Note that exporting a variable from a Jupyter *terminal* does **not** reach the
notebook. Both are children of the Jupyter server, which captured its
environment when it started, so a kernel restart re-inherits the same
environment and still cannot see it. `get_workspace_bucket` therefore does not
rely on the variable alone -- see its docstring.

None of these functions is required to use ``aoutools`` -- ``calculate_prs``
and friends take a VDS you have already read. They exist to keep the
boilerplate short and correct.
"""

import json
import logging
import os
import subprocess
import warnings

import hail as hl

logger = logging.getLogger(__name__)

# The Verily Workbench CLI. It is the authoritative source for what this
# workspace's resources actually are, which the environment and the mount table
# are only proxies for. Absent off the Workbench.
_WB_CLI = "wb"
_WB_TIMEOUT_SECONDS = 30

# The workspace bucket's resource `id` contains "workspace-bucket" -- but so
# does "temporary-workspace-bucket", whose contents are garbage-collected. The
# temp bucket must therefore be excluded *before* the substring match, not
# after. Verily's own setup notebook only gets this right because it happens to
# test `temporary` first.
_WORKSPACE_BUCKET_ID = "workspace-bucket"
_TEMP_BUCKET_MARKER = "temporary"

# Where gcsfuse mounts show up in /proc/mounts. The *device* field of a gcsfuse
# entry is the bucket name. This is a **fallback** for when the CLI is
# unavailable: a workspace can have several buckets mounted (a cloned tutorial
# workspace brings the notebooks bucket it was cloned from), so the workspace
# bucket is identified by its mountpoint, never by position in the table.
_MOUNTS_FILE = "/proc/mounts"
_GCSFUSE_TYPES = ("fuse.gcsfuse", "gcsfuse")
_WORKSPACE_MOUNT_NAME = "workspace-bucket"


def _wb_json(*args: str) -> object | None:
    """
    Runs a `wb` CLI command and parses its JSON, or returns None.

    Returns None -- rather than raising -- whenever the CLI is missing, times
    out, fails, or emits something that is not JSON. Every caller has another
    source to fall back on, and `wb` does not exist off the Workbench at all.
    """
    try:
        result = subprocess.run(
            [_WB_CLI, *args, "--format=json"],
            capture_output=True,
            text=True,
            check=True,
            timeout=_WB_TIMEOUT_SECONDS,
        )
        return json.loads(result.stdout)
    except (
        OSError,
        subprocess.SubprocessError,
        json.JSONDecodeError,
    ) as e:
        logger.debug("`wb %s` did not yield JSON: %s", " ".join(args), e)
        return None


def _bucket_from_wb_cli() -> str | None:
    """
    Asks the Verily Workbench CLI for this workspace's bucket, or None.

    This is the **authoritative** source: it reports the workspace's declared
    resources, rather than inferring them from an environment variable that may
    not be exported or a mount table that may hold several buckets.

    Selects the `GCS_BUCKET` resource whose `id` names the workspace bucket,
    explicitly skipping the *temporary* workspace bucket -- whose id also
    contains "workspace-bucket", and whose contents are garbage-collected.
    """
    resources = _wb_json("resource", "list")
    if not isinstance(resources, list):
        return None

    for resource in resources:
        if not isinstance(resource, dict):
            continue
        if resource.get("resourceType") != "GCS_BUCKET":
            continue

        resource_id = str(resource.get("id", ""))
        if _TEMP_BUCKET_MARKER in resource_id:
            continue  # the temp bucket is GC'd; never write results there
        if _WORKSPACE_BUCKET_ID not in resource_id:
            continue  # e.g. `aou-tutorial-notebooks`, a cloned workspace's

        bucket_name = resource.get("bucketName")
        if bucket_name:
            return str(bucket_name)
    return None


def _project_from_wb_cli() -> str | None:
    """Asks the Workbench CLI for this workspace's Google project, or None."""
    workspace = _wb_json("workspace", "describe")
    if isinstance(workspace, dict):
        project = workspace.get("googleProjectId")
        if project:
            return str(project)
    return None


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


def _gcsfuse_mounts() -> list[tuple[str, str]]:
    """
    Every gcsfuse mount, as a list of `(bucket_name, mountpoint)`.

    In `/proc/mounts` the *device* field of a gcsfuse entry is the bucket name.
    Returns an empty list on any platform without `/proc/mounts` (macOS,
    Windows), or if nothing is fuse-mounted.
    """
    try:
        with open(_MOUNTS_FILE) as fh:
            lines = fh.readlines()
    except OSError:
        return []

    mounts = []
    for line in lines:
        fields = line.split()
        if len(fields) < 3:
            continue
        device, mountpoint, fstype = fields[0], fields[1], fields[2]
        if fstype in _GCSFUSE_TYPES and device:
            mounts.append((device, mountpoint))
    return mounts


def _bucket_from_gcsfuse_mount() -> str | None:
    """
    Recovers the workspace bucket from the gcsfuse mount table, or None.

    **A workspace can have more than one bucket fuse-mounted.** A cloned
    tutorial workspace, for instance, mounts the notebooks bucket it was cloned
    from alongside the workspace's own bucket. So this cannot simply take the
    first gcsfuse entry it finds: doing that once selected
    `cloned-aou-tutorial-notebooks-...` in place of `workspace-bucket-...`, and
    staged a user's files into the wrong bucket without a word.

    The workspace bucket is the one mounted at `.../workspace-bucket`, so it is
    identified by its **mountpoint**, not by position in the table.

    Returns None if the answer is not unambiguous -- if no gcsfuse mount matches
    the expected mountpoint and there is more than one candidate, the caller
    raises rather than picking one. Writing analysis output to the wrong bucket
    is a silent error, and a guess is not worth making.
    """
    mounts = _gcsfuse_mounts()
    if not mounts:
        return None

    # The workspace bucket lives at `~/workspace/workspace-bucket`.
    for bucket, mountpoint in mounts:
        if os.path.basename(mountpoint.rstrip("/")) == _WORKSPACE_MOUNT_NAME:
            return bucket

    # No mountpoint says which one it is. If there is exactly one bucket
    # mounted, it is not a guess. Otherwise, refuse.
    if len(mounts) == 1:
        return mounts[0][0]

    return None


def get_workspace_bucket(bucket: str | None = None) -> str:
    """
    Resolves the GCS path of the workspace bucket, e.g. `gs://my-bucket`.

    Checks four sources, in order, and returns the first that yields an answer:

    1. The `bucket` argument, if given.
    2. The `WORKSPACE_BUCKET` environment variable.
    3. The **Verily Workbench CLI** (`wb resource list`), which reports the
       workspace's declared resources. This is the authoritative source, and it
       is what Verily's own setup notebook uses.
    4. The **gcsfuse mount table**, as a fallback if the CLI is unavailable.

    Why not just tell users to export the variable: **exporting it from a
    Jupyter terminal does not work.** The terminal and the kernel are sibling
    children of the Jupyter server, which captured its environment at startup,
    so the kernel never sees a variable exported later in a terminal -- and
    restarting the kernel re-inherits that same, unchanged environment.

    To set it from inside a notebook (which *does* work, since it mutates the
    kernel's own environment)::

        import os
        os.environ["WORKSPACE_BUCKET"] = "gs://your-bucket"

    To make that persist across kernel restarts, drop the same two lines into
    `~/.ipython/profile_default/startup/00-aou-env.py`, which every kernel runs
    at startup.

    Notes
    -----
    A workspace commonly has **more than one** bucket, and picking the wrong one
    does not fail -- the results simply are not where you think they are. Two in
    particular are excluded here:

    * the **temporary** workspace bucket, whose resource id also contains
      "workspace-bucket" and whose contents are garbage-collected;
    * a cloned workspace's notebooks bucket (e.g. `aou-tutorial-notebooks`),
      which is fuse-mounted alongside your own.

    Parameters
    ----------
    bucket : str, optional
        An explicit bucket, with or without the `gs://` prefix. Overrides every
        other source.

    Returns
    -------
    str
        The bucket as a GCS path, with no trailing slash (e.g. `gs://my-bucket`).

    Raises
    ------
    OSError
        If the bucket cannot be determined, or if several are mounted and none
        can be identified as the workspace bucket. It will not guess.
    """
    if not bucket:
        bucket = os.getenv("WORKSPACE_BUCKET")

    if not bucket:
        from_cli = _bucket_from_wb_cli()
        if from_cli:
            bucket = from_cli
            logger.info(
                "Resolved the workspace bucket from the Workbench CLI: gs://%s",
                from_cli,
            )

    if not bucket:
        discovered = _bucket_from_gcsfuse_mount()
        if discovered:
            bucket = discovered
            warnings.warn(
                "The 'WORKSPACE_BUCKET' environment variable is not set, so "
                "aoutools recovered the workspace bucket from the gcsfuse "
                f"mount table:\n  gs://{discovered}\n"
                "Check that this is the bucket you expect. To set it "
                "explicitly, run this in a notebook cell (exporting it from a "
                "terminal will NOT reach the kernel):\n"
                '  import os; os.environ["WORKSPACE_BUCKET"] = "gs://..."',
                UserWarning,
                stacklevel=2,
            )
        else:
            # Several buckets are mounted and none of them sits at the
            # workspace-bucket mountpoint. Do NOT pick one: output written to
            # the wrong bucket is a silent error, and the candidates here are
            # things like a cloned tutorial's notebooks bucket, which a user
            # would never intend to write results into.
            candidates = _gcsfuse_mounts()
            if len(candidates) > 1:
                listed = "\n".join(
                    f"  gs://{b}   (mounted at {m})" for b, m in candidates
                )
                raise OSError(
                    "Several buckets are fuse-mounted and none of them is at "
                    "the expected workspace-bucket mountpoint, so aoutools "
                    "cannot tell which one is yours:\n\n"
                    f"{listed}\n\n"
                    "It will not guess: writing results into the wrong bucket "
                    "would fail silently. Say which one, in a notebook cell:\n"
                    "  import os\n"
                    '  os.environ["WORKSPACE_BUCKET"] = "gs://your-bucket"'
                )

    if not bucket:
        raise OSError(
            "Could not determine the workspace bucket. The "
            "'WORKSPACE_BUCKET' environment variable is not set, and no "
            "gcsfuse mount was found.\n\n"
            "Current All of Us (Verily) images do not export this variable. "
            "Note that exporting it from a Jupyter *terminal* does not reach "
            "the notebook kernel -- they are sibling processes, so the kernel "
            "never sees it, even after a restart.\n\n"
            "Set it from inside a notebook cell instead:\n"
            "  import os\n"
            '  os.environ["WORKSPACE_BUCKET"] = "gs://your-bucket"\n\n'
            "To persist it across kernel restarts, put those lines in\n"
            "  ~/.ipython/profile_default/startup/00-aou-env.py"
        )

    bucket = bucket.rstrip("/")
    if not bucket.startswith("gs://"):
        bucket = f"gs://{bucket}"
    return bucket


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


def get_google_project(project: str | None = None) -> str | None:
    """
    Resolves the Google project to bill requester-pays reads to.

    Checks, in order: the `project` argument, `$GOOGLE_PROJECT`,
    `$GOOGLE_CLOUD_PROJECT`, and finally the Workbench CLI
    (`wb workspace describe`).

    Two variable names are checked because the platforms disagree. The older
    All of Us Workbench exported `GOOGLE_PROJECT`; Verily's setup notebook sets
    `GOOGLE_CLOUD_PROJECT` and leaves `GOOGLE_PROJECT` empty. Reading only the
    first meant a clean Verily image silently got no billing project at all,
    and every read of the requester-pays VDS bucket failed.

    Parameters
    ----------
    project : str, optional
        An explicit project, which overrides every other source.

    Returns
    -------
    str or None
        The project id, or None if no source knows one (e.g. off the Workbench).
    """
    if project:
        return project

    for var in ("GOOGLE_PROJECT", "GOOGLE_CLOUD_PROJECT"):
        value = os.getenv(var)
        if value:
            return value

    from_cli = _project_from_wb_cli()
    if from_cli:
        logger.info(
            "Resolved the Google project from the Workbench CLI: %s", from_cli
        )
    return from_cli


def init_hail(reference: str = "GRCh38", **kwargs) -> None:
    """
    Initializes Hail for the All of Us Researcher Workbench.

    Wraps `hl.init` with the two settings the Workbench needs:

    * **Requester-pays billing.** The bucket holding the VDS charges the reader,
      so Hail must be given a project to bill; without one, reading the VDS
      fails. The project is resolved by `get_google_project`, which checks
      `GOOGLE_PROJECT`, then `GOOGLE_CLOUD_PROJECT`, then the Workbench CLI --
      the older Workbench exported the first, the current (Verily) platform uses
      the second, and a clean image may set neither.
    * **The reference genome.** `hl.init(default_reference=...)` is deprecated;
      the reference is now set by calling `hl.default_reference` *with* an
      argument, after initialization. Doing it here keeps the deprecation
      warning out of your notebook.

    Calling this off the Workbench, where no project can be found, is fine: the
    requester-pays setting is simply omitted, with a warning.

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

    google_project = get_google_project()
    if google_project:
        kwargs.setdefault("gcs_requester_pays_configuration", google_project)
    else:
        warnings.warn(
            "No Google project could be found (GOOGLE_PROJECT, "
            "GOOGLE_CLOUD_PROJECT, and the Workbench CLI were all tried), so "
            "Hail is being initialized WITHOUT a requester-pays billing "
            "project. Reading the All of Us VDS will fail with a permissions "
            "error. Pass one explicitly if you are on the Workbench:\n"
            '  init_hail(gcs_requester_pays_configuration="your-project")',
            UserWarning,
            stacklevel=2,
        )

    hl.init(**kwargs)
    # Must come after hl.init: passing default_reference to hl.init is
    # deprecated, and calling hl.default_reference() with no argument is the
    # getter. With an argument it is the setter.
    hl.default_reference(reference)
