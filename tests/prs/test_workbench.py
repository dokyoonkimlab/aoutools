"""Tests for the Workbench environment helpers.

These pin down three things the Workbench changed under us, each of which every
notebook would otherwise get wrong in its own way:

* `hl.init(default_reference=...)` is deprecated. The reference genome must be
  set *after* init, by calling `hl.default_reference` **with** an argument.
* `WGS_VDS_PATH` is no longer exported, so the VDS path has to fall back to a
  pinned default -- and that fallback must be loud, because it names a specific
  All of Us data release that this package cannot check against the workspace.
* `WORKSPACE_BUCKET` is no longer exported *either*, and it cannot simply be
  exported by the user: a variable set in a Jupyter terminal never reaches the
  kernel, because the two are sibling children of the Jupyter server. So the
  bucket is recovered from the gcsfuse mount table instead.
"""

import warnings
from unittest.mock import mock_open, patch

import pytest

from aoutools._workbench import (
    DEFAULT_VDS_PATH,
    get_vds_path,
    get_workspace_bucket,
    init_hail,
)

# A realistic /proc/mounts, with the workspace bucket fuse-mounted.
MOUNTS = (
    "proc /proc proc rw,nosuid,nodev,noexec,relatime 0 0\n"
    "/dev/sda1 / ext4 rw,relatime 0 0\n"
    "workspace-bucket-wb-swift-orange-3552 "
    "/home/dataproc/workspace/workspace-bucket fuse.gcsfuse "
    "rw,nosuid,nodev,relatime,user_id=1000 0 0\n"
)
MOUNTS_NO_FUSE = (
    "proc /proc proc rw,nosuid,nodev,noexec,relatime 0 0\n"
    "/dev/sda1 / ext4 rw,relatime 0 0\n"
)

# The layout that broke it. A cloned tutorial workspace mounts the notebooks
# bucket it was cloned from ALONGSIDE the workspace's own bucket -- and the
# kernel happened to list the cloned one first. Taking the first gcsfuse entry
# selected `cloned-aou-tutorial-notebooks-...` and staged a user's files into
# the wrong bucket without a word.
MOUNTS_MULTI = (
    "proc /proc proc rw,nosuid,nodev,noexec,relatime 0 0\n"
    "cloned-aou-tutorial-notebooks-wb-swift-orange-3552 "
    "/home/dataproc/workspace/cloned-aou-tutorial-notebooks fuse.gcsfuse "
    "rw,nosuid,nodev,relatime 0 0\n"
    "workspace-bucket-wb-swift-orange-3552 "
    "/home/dataproc/workspace/workspace-bucket fuse.gcsfuse "
    "rw,nosuid,nodev,relatime 0 0\n"
)

# Several buckets mounted, none of them at the workspace-bucket mountpoint.
MOUNTS_AMBIGUOUS = (
    "proc /proc proc rw,nosuid,nodev,noexec,relatime 0 0\n"
    "bucket-one /home/dataproc/workspace/one fuse.gcsfuse rw 0 0\n"
    "bucket-two /home/dataproc/workspace/two fuse.gcsfuse rw 0 0\n"
)


class TestGetWorkspaceBucket:
    """Resolution: explicit argument, then env var, then the gcsfuse mount."""

    def test_explicit_bucket_wins(self, monkeypatch):
        monkeypatch.setenv("WORKSPACE_BUCKET", "gs://from-env")

        assert get_workspace_bucket("gs://explicit") == "gs://explicit"

    def test_adds_the_gs_prefix_and_strips_a_trailing_slash(self):
        assert get_workspace_bucket("my-bucket/") == "gs://my-bucket"

    def test_env_var_is_used_when_set(self, monkeypatch):
        monkeypatch.setenv("WORKSPACE_BUCKET", "gs://from-env")

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert get_workspace_bucket() == "gs://from-env"

    def test_recovers_the_bucket_from_the_gcsfuse_mount(self, monkeypatch):
        """The case that matters on the current (Verily) platform.

        `WORKSPACE_BUCKET` is not exported, and the user cannot fix that from a
        terminal -- the kernel is a sibling process and never sees it. But the
        bucket is fuse-mounted into the workspace regardless, and gcsfuse puts
        the bucket name in the *device* field of the mount table.
        """
        monkeypatch.delenv("WORKSPACE_BUCKET", raising=False)

        with patch("builtins.open", mock_open(read_data=MOUNTS)):
            with pytest.warns(UserWarning, match="gcsfuse"):
                bucket = get_workspace_bucket()

        assert bucket == "gs://workspace-bucket-wb-swift-orange-3552"

    def test_picks_the_workspace_bucket_not_whichever_is_mounted_first(
        self, monkeypatch
    ):
        """The regression. A cloned tutorial workspace has TWO buckets mounted.

        Taking the first gcsfuse entry in the table selected the cloned
        tutorial's notebooks bucket and quietly staged the user's weights file
        into it. The workspace bucket is identified by its *mountpoint*
        (`.../workspace-bucket`), not by its position in the mount table.
        """
        monkeypatch.delenv("WORKSPACE_BUCKET", raising=False)

        with patch("builtins.open", mock_open(read_data=MOUNTS_MULTI)):
            with pytest.warns(UserWarning):
                bucket = get_workspace_bucket()

        assert bucket == "gs://workspace-bucket-wb-swift-orange-3552"
        assert "cloned" not in bucket, (
            "the cloned tutorial's notebooks bucket is not the workspace bucket"
        )

    def test_refuses_to_guess_between_several_buckets(self, monkeypatch):
        """When nothing identifies the workspace bucket, do not pick one.

        Output written to the wrong bucket fails silently -- there is no error,
        the results simply are not where the user thinks they are. An error
        naming the candidates is strictly better than a coin flip.
        """
        monkeypatch.delenv("WORKSPACE_BUCKET", raising=False)

        with patch("builtins.open", mock_open(read_data=MOUNTS_AMBIGUOUS)):
            with pytest.raises(OSError) as excinfo:
                get_workspace_bucket()

        message = str(excinfo.value)
        assert "bucket-one" in message and "bucket-two" in message, (
            "the error must name the candidates it found"
        )
        assert "will not guess" in message

    def test_a_single_mounted_bucket_is_not_a_guess(self, monkeypatch):
        """One bucket, non-standard mountpoint: unambiguous, so use it."""
        monkeypatch.delenv("WORKSPACE_BUCKET", raising=False)
        mounts = (
            "only-bucket /home/dataproc/somewhere-else fuse.gcsfuse rw 0 0\n"
        )

        with patch("builtins.open", mock_open(read_data=mounts)):
            with pytest.warns(UserWarning):
                assert get_workspace_bucket() == "gs://only-bucket"

    def test_raises_with_instructions_when_nothing_works(self, monkeypatch):
        """The error has to say how to fix it, because the obvious fix (export
        it from a terminal) is the one thing that does NOT work."""
        monkeypatch.delenv("WORKSPACE_BUCKET", raising=False)

        with patch("builtins.open", mock_open(read_data=MOUNTS_NO_FUSE)):
            with pytest.raises(OSError) as excinfo:
                get_workspace_bucket()

        message = str(excinfo.value)
        assert "os.environ" in message, "must show the fix that works"
        assert "terminal" in message, "must warn off the fix that does not"

    def test_missing_proc_mounts_is_not_a_crash(self, monkeypatch):
        """No /proc/mounts on macOS. Must degrade to the clear error, not an
        OSError from the open() itself."""
        monkeypatch.delenv("WORKSPACE_BUCKET", raising=False)

        with patch("builtins.open", side_effect=FileNotFoundError):
            with pytest.raises(OSError, match="Could not determine"):
                get_workspace_bucket()


class TestGetVdsPath:
    """Resolution order: explicit argument, then env var, then the default."""

    def test_explicit_path_wins_over_everything(self, monkeypatch):
        monkeypatch.setenv("WGS_VDS_PATH", "gs://from-env/hail.vds")

        assert get_vds_path("gs://explicit/hail.vds") == (
            "gs://explicit/hail.vds"
        )

    def test_env_var_is_used_when_no_path_is_given(self, monkeypatch):
        monkeypatch.setenv("WGS_VDS_PATH", "gs://from-env/hail.vds")

        # No warning: the environment is authoritative, so there is nothing to
        # caution the user about.
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert get_vds_path() == "gs://from-env/hail.vds"

    def test_falls_back_to_the_default_and_says_so(self, monkeypatch):
        """The fallback must warn.

        It names a specific data release (v9 today). If the workspace is
        registered against a different one, the genotypes would silently
        disagree with the phenotypes -- a wrong analysis, not an error. The
        warning is the only thing standing between the user and that.
        """
        monkeypatch.delenv("WGS_VDS_PATH", raising=False)

        with pytest.warns(UserWarning, match="WGS_VDS_PATH"):
            assert get_vds_path() == DEFAULT_VDS_PATH

    def test_empty_env_var_is_treated_as_unset(self, monkeypatch):
        monkeypatch.setenv("WGS_VDS_PATH", "")

        with pytest.warns(UserWarning):
            assert get_vds_path() == DEFAULT_VDS_PATH


class TestInitHail:
    """`init_hail` wraps `hl.init` with the two Workbench-specific settings."""

    def test_sets_the_reference_after_init_not_through_it(self, monkeypatch):
        """The whole reason this helper exists.

        `hl.init(default_reference=...)` still works but is deprecated, and
        prints a warning into every notebook. The replacement is to call
        `hl.default_reference` with an argument once Hail is up.
        """
        monkeypatch.setenv("GOOGLE_PROJECT", "my-project")

        with patch("aoutools._workbench.hl") as mock_hl:
            init_hail()

        _, init_kwargs = mock_hl.init.call_args
        assert "default_reference" not in init_kwargs, (
            "passing default_reference to hl.init is deprecated"
        )
        mock_hl.default_reference.assert_called_once_with("GRCh38")

    def test_bills_reads_to_the_workspace_project(self, monkeypatch):
        """Without this, reading the requester-pays VDS bucket fails."""
        monkeypatch.setenv("GOOGLE_PROJECT", "my-project")

        with patch("aoutools._workbench.hl") as mock_hl:
            init_hail()

        _, init_kwargs = mock_hl.init.call_args
        assert init_kwargs["gcs_requester_pays_configuration"] == "my-project"
        assert init_kwargs["idempotent"] is True

    def test_works_off_the_workbench(self, monkeypatch):
        """No GOOGLE_PROJECT means no requester-pays argument, not a KeyError.

        Reading the All of Us VDS will fail in this state, but everything else
        (local hail, the test suites) works, so the helper must not hard-fail.
        """
        monkeypatch.delenv("GOOGLE_PROJECT", raising=False)

        with patch("aoutools._workbench.hl") as mock_hl:
            init_hail()

        _, init_kwargs = mock_hl.init.call_args
        assert "gcs_requester_pays_configuration" not in init_kwargs
        mock_hl.default_reference.assert_called_once_with("GRCh38")

    def test_caller_can_override_the_defaults(self, monkeypatch):
        monkeypatch.setenv("GOOGLE_PROJECT", "my-project")

        with patch("aoutools._workbench.hl") as mock_hl:
            init_hail(
                reference="GRCh37",
                idempotent=False,
                gcs_requester_pays_configuration="other-project",
                master="local[1]",
            )

        _, init_kwargs = mock_hl.init.call_args
        assert init_kwargs["idempotent"] is False
        assert init_kwargs["gcs_requester_pays_configuration"] == (
            "other-project"
        )
        assert init_kwargs["master"] == "local[1]"
        mock_hl.default_reference.assert_called_once_with("GRCh37")
