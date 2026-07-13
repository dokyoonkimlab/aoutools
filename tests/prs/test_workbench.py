"""Tests for the Workbench environment helpers.

These pin down two things the Workbench changed under us, and which every
notebook would otherwise get wrong in its own way:

* `hl.init(default_reference=...)` is deprecated. The reference genome must be
  set *after* init, by calling `hl.default_reference` **with** an argument.
* `WGS_VDS_PATH` is no longer exported, so the VDS path has to fall back to a
  pinned default -- and that fallback must be loud, because it names a specific
  All of Us data release that this package cannot check against the workspace.
"""

import warnings
from unittest.mock import patch

import pytest

from aoutools._workbench import DEFAULT_VDS_PATH, get_vds_path, init_hail


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
