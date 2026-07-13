"""Tests for the PGS Catalog downloader's environment handling.

`_downloader.py` shells out to a `pgscatalog-download` CLI that it installs
into an isolated venv on first use. Both steps can fail, and both used to fail
*silently* -- which is the entire subject of this file:

* `_run` logged the subprocess's stderr through `logger.error`. `aoutools`
  installs a `NullHandler` and a notebook configures no logging, so the reason
  for a failed download went nowhere and the caller saw a bare exit code.
* `_ensure_pgscatalog_download` cached and returned the CLI path without ever
  checking the file existed, so a half-finished install was cached as a success.
"""

import subprocess
from pathlib import Path

import pytest

from aoutools.prs import _downloader
from aoutools.prs._downloader import _ensure_pgscatalog_download, _run


@pytest.fixture(autouse=True)
def clear_cli_cache():
    """The CLI path is cached in a module global; reset it between tests."""
    _downloader._cached_cli_path = None
    yield
    _downloader._cached_cli_path = None


class TestRun:
    def test_failure_carries_the_subprocess_stderr(self, mocker):
        """The error text must reach the *exception*, not just the logger.

        Logging is a NullHandler by default, so a caller who only sees the
        exception used to get 'returned non-zero exit status 1' and nothing
        else -- no way to tell a missing directory from a network failure.
        """
        mocker.patch(
            "aoutools.prs._downloader.subprocess.run",
            side_effect=subprocess.CalledProcessError(
                returncode=1,
                cmd=["pgscatalog-download", "-o", "/nope"],
                output="",
                stderr="FileNotFoundError: --outdir nope doesn't exist",
            ),
        )

        with pytest.raises(RuntimeError) as excinfo:
            _run(["pgscatalog-download", "-o", "/nope"])

        message = str(excinfo.value)
        assert "--outdir nope doesn't exist" in message
        assert "exit code 1" in message

    def test_failure_keeps_the_original_exception_as_the_cause(self, mocker):
        """Chained, so the CalledProcessError is still reachable."""
        original = subprocess.CalledProcessError(
            returncode=2, cmd=["x"], output="", stderr="boom"
        )
        mocker.patch(
            "aoutools.prs._downloader.subprocess.run", side_effect=original
        )

        with pytest.raises(RuntimeError) as excinfo:
            _run(["x"])

        assert excinfo.value.__cause__ is original

    def test_success_is_quiet(self, mocker):
        mocker.patch(
            "aoutools.prs._downloader.subprocess.run",
            return_value=subprocess.CompletedProcess(["x"], 0, "ok", ""),
        )

        assert _run(["x"]) is None


class TestEnsurePgscatalogDownload:
    def test_raises_when_the_cli_is_missing_after_install(
        self, mocker, tmp_path
    ):
        """A half-built venv must fail here, loudly, not later.

        Nothing verifies that `pip install` produced a working CLI. Without
        this check the missing path is cached as a success, and the user meets
        it as an inscrutable failure from a command that isn't there.
        """
        mocker.patch.object(_downloader, "PGS_ENV_DIR", tmp_path / "env")
        bin_dir = tmp_path / "env" / "bin"
        bin_dir.mkdir(parents=True)
        (bin_dir / "pip").touch()  # venv looks built...
        # ...but pgscatalog-download was never produced.
        mocker.patch.object(
            _downloader, "_get_pgscatalog_version", return_value=None
        )
        mocker.patch.object(_downloader, "_run")  # pretend pip "succeeded"

        with pytest.raises(RuntimeError, match="is not there"):
            _ensure_pgscatalog_download()

        assert _downloader._cached_cli_path is None, (
            "a failed lookup must not be cached as success"
        )

    def test_returns_and_caches_the_cli_when_present(self, mocker, tmp_path):
        mocker.patch.object(_downloader, "PGS_ENV_DIR", tmp_path / "env")
        bin_dir = tmp_path / "env" / "bin"
        bin_dir.mkdir(parents=True)
        (bin_dir / "pip").touch()
        (bin_dir / "pgscatalog-download").touch()
        mocker.patch.object(
            _downloader, "_get_pgscatalog_version", return_value="1.0.1"
        )
        mock_run = mocker.patch.object(_downloader, "_run")

        cli = _ensure_pgscatalog_download(min_version="1.0.1")

        assert cli == bin_dir / "pgscatalog-download"
        assert isinstance(cli, Path)
        mock_run.assert_not_called()  # already at min_version: no reinstall
        assert _downloader._cached_cli_path == cli
