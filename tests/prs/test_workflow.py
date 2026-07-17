"""Unit tests for the `_workflow.py` end-to-end PGS workflow.

Hail and the download/read/score steps are mocked; these tests pin the control
flow and the operator-facing WARNING emitted when one scoring file cannot be
read but the rest of the batch can (see the logging policy in AGENTS.md).
"""

import logging
from pathlib import Path
from unittest.mock import MagicMock

from aoutools.prs._workflow import calculate_pgs


def _write_downloads(names):
    """A `download_pgs` side effect that drops empty files into `outdir`."""

    def _side_effect(*, outdir, **kwargs):
        for name in names:
            Path(outdir, name).touch()

    return _side_effect


def test_unreadable_scoring_file_is_skipped_with_warning(mocker, caplog):
    """One bad file warns and is skipped; the batch runs on the survivor.

    The failure mode this guards against is the opposite: a single unreadable
    score taking down the whole run silently. It must WARN, name the score, and
    keep going.
    """
    mocker.patch(
        "aoutools.prs._workflow.download_pgs",
        side_effect=_write_downloads(["PGS000001.txt.gz", "PGS000002.txt.gz"]),
    )

    def _fake_read(**kwargs):
        if "PGS000001" in kwargs["file_path"]:
            raise ValueError("scoring file has no positional information")
        return MagicMock(name="weights_table")

    mocker.patch(
        "aoutools.prs._workflow.read_prs_weights", side_effect=_fake_read
    )
    batch = mocker.patch(
        "aoutools.prs._workflow.calculate_prs_batch",
        return_value="gs://bucket/out.csv",
    )

    with caplog.at_level(logging.WARNING, logger="aoutools.prs._workflow"):
        result = calculate_pgs(
            vds=MagicMock(),
            output_path="gs://bucket/out.csv",
            pgs=["PGS000001", "PGS000002"],
        )

    assert result == "gs://bucket/out.csv"
    assert "Could not read scoring file for 'PGS000001'" in caplog.text
    assert any(r.levelno == logging.WARNING for r in caplog.records)

    # Only the readable score reaches the batch step.
    _, kwargs = batch.call_args
    assert set(kwargs["weights_tables_map"]) == {"PGS000002"}


def test_all_files_unreadable_returns_none(mocker, caplog):
    """If nothing can be read, warn and return None -- never call the scorer."""
    mocker.patch(
        "aoutools.prs._workflow.download_pgs",
        side_effect=_write_downloads(["PGS000001.txt.gz"]),
    )
    mocker.patch(
        "aoutools.prs._workflow.read_prs_weights",
        side_effect=ValueError("bad file"),
    )
    batch = mocker.patch("aoutools.prs._workflow.calculate_prs_batch")

    with caplog.at_level(logging.WARNING, logger="aoutools.prs._workflow"):
        result = calculate_pgs(
            vds=MagicMock(),
            output_path="gs://bucket/out.csv",
            pgs="PGS000001",
        )

    assert result is None
    batch.assert_not_called()
    assert "No valid scoring files could be read" in caplog.text
