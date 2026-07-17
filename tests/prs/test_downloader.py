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
from aoutools.prs._downloader import (
    _ensure_pgscatalog_download,
    _run,
    _run_pgscatalog_download,
)


@pytest.fixture(autouse=True)
def clear_cli_cache():
    """The CLI path is cached in a module global; reset it between tests."""
    _downloader._cached_cli_path = None
    yield
    _downloader._cached_cli_path = None


class TestSubprocessEnv:
    """The isolated venv only isolates imports; the child processes must also
    run with a scrubbed environment, or the base image's pinned tenacity
    (8.2.3, held by dsub) defeats the >=9.0.0 pgscatalog.core needs -- at
    install time via PIP_CONSTRAINT, or at runtime via PYTHONPATH."""

    def test_strips_vars_that_defeat_the_venv(self, monkeypatch):
        monkeypatch.setenv("PIP_CONSTRAINT", "/etc/pip/constraints.txt")
        monkeypatch.setenv(
            "PYTHONPATH", "/opt/conda/lib/python3.10/site-packages"
        )
        monkeypatch.setenv("PYTHONHOME", "/opt/conda")
        monkeypatch.setenv("PIP_USER", "1")

        env = _downloader._subprocess_env()

        assert "PIP_CONSTRAINT" not in env
        assert "PYTHONPATH" not in env
        assert "PYTHONHOME" not in env
        assert "PIP_USER" not in env
        assert env["PYTHONNOUSERSITE"] == "1"

    def test_leaves_the_package_index_alone(self, monkeypatch):
        """A Workbench mirror is load-bearing for the download itself."""
        monkeypatch.setenv("PIP_INDEX_URL", "https://mirror.internal/simple")

        env = _downloader._subprocess_env()

        assert env["PIP_INDEX_URL"] == "https://mirror.internal/simple"

    def test_run_applies_the_scrubbed_env_by_default(self, mocker, monkeypatch):
        monkeypatch.setenv("PIP_CONSTRAINT", "/etc/pip/constraints.txt")
        run = mocker.patch(
            "aoutools.prs._downloader.subprocess.run",
            return_value=subprocess.CompletedProcess(["x"], 0, "ok", ""),
        )

        _run(["x"])

        passed_env = run.call_args.kwargs["env"]
        assert "PIP_CONSTRAINT" not in passed_env
        assert passed_env["PYTHONNOUSERSITE"] == "1"


class TestGetPgscatalogVersion:
    """The version probe runs `pip show` before install, when the package is
    legitimately absent. `pip show` then prints 'Package(s) not found' to
    stderr -- expected, but it reads like a fault on every first run, so the
    probe swallows stderr while still returning None on the non-zero exit."""

    def test_absent_package_returns_none_without_leaking_stderr(self, mocker):
        check_output = mocker.patch(
            "aoutools.prs._downloader.subprocess.check_output",
            side_effect=subprocess.CalledProcessError(
                returncode=1, cmd=["pip", "show", "pgscatalog.core"]
            ),
        )

        result = _downloader._get_pgscatalog_version(Path("/env/bin"))

        assert result is None
        assert check_output.call_args.kwargs["stderr"] == subprocess.DEVNULL, (
            "pip's 'not found' stderr must not reach the notebook"
        )

    def test_parses_the_installed_version(self, mocker):
        mocker.patch(
            "aoutools.prs._downloader.subprocess.check_output",
            return_value="Name: pgscatalog.core\nVersion: 1.0.1\n",
        )

        assert _downloader._get_pgscatalog_version(Path("/env/bin")) == "1.0.1"


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


class TestRunPgscatalogDownload:
    def test_an_existing_file_gets_an_actionable_error(self, mocker, tmp_path):
        """A pre-existing scoring file must not surface as a raw traceback.

        Without `-w`, pgscatalog-download raises FileExistsError on the first
        file already present -- and the downloads are concurrent, so that takes
        down the *whole batch*, including scores that would have worked.
        Re-running a notebook cell is the ordinary way to hit this, and the
        underlying error arrives buried inside tenacity and a thread pool with
        no hint of the fix.
        """
        mocker.patch.object(
            _downloader,
            "_ensure_pgscatalog_download",
            return_value=tmp_path / "pgscatalog-download",
        )
        mocker.patch.object(
            _downloader,
            "_run",
            side_effect=RuntimeError(
                "Command failed with exit code 1:\n"
                "FileExistsError: /out/PGS000747_hmPOS_GRCh38.txt.gz "
                "already exists"
            ),
        )

        with pytest.raises(FileExistsError) as excinfo:
            _run_pgscatalog_download("/out", "-i", "PGS000746", "PGS000747")

        message = str(excinfo.value)
        assert "overwrite_existing_file=True" in message, (
            "the error must name the flag that fixes it"
        )
        assert "entire batch" in message

    def test_other_failures_are_not_disguised_as_file_conflicts(
        self, mocker, tmp_path
    ):
        """Only the already-exists case gets rewritten. A network failure must
        still arrive as itself."""
        mocker.patch.object(
            _downloader,
            "_ensure_pgscatalog_download",
            return_value=tmp_path / "pgscatalog-download",
        )
        mocker.patch.object(
            _downloader,
            "_run",
            side_effect=RuntimeError("Connection reset by peer"),
        )

        with pytest.raises(RuntimeError, match="Connection reset"):
            _run_pgscatalog_download("/out", "-i", "PGS000746")


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
