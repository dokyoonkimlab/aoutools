"""Chunked scoring, against a real hail VDS.

The chunking strategy is the reason the library is affordable on the All of Us
VDS: the weights are cut into `chunk_size` pieces, each piece builds 1bp
intervals so `hl.vds.filter_intervals` reads only the loci it needs, and the
per-chunk scores are summed. That summation is only correct if the chunks
partition the weights exactly -- an overlap double-counts a variant, a gap drops
one, and either way the score is still a clean float.

Nothing here is mocked except `hfs.open`, which would otherwise require GCS.
"""

import pandas as pd
import pytest

from aoutools.prs import _calculator
from aoutools.prs._calculator import _aggregate_and_export, _process_chunks
from aoutools.prs._calculator_utils import _prepare_weights_for_chunking
from aoutools.prs._config import PRSConfig

pytestmark = pytest.mark.integration


def totals(vds, raw_weights, chunk_size, **config_kwargs):
    """Scores the VDS in chunks and returns {sample: summed prs}."""
    config = PRSConfig(chunk_size=chunk_size, **config_kwargs)
    weights, n_chunks = _prepare_weights_for_chunking(
        weights_table=raw_weights, config=config, validate_table=True
    )
    partial_dfs = _process_chunks(
        full_weights_table=weights,
        n_chunks=n_chunks,
        vds=vds,
        config=config,
    )
    combined = pd.concat(partial_dfs, ignore_index=True)
    summed = combined.groupby("person_id").sum()
    return n_chunks, summed["prs"].to_dict()


@pytest.mark.parametrize("split_multi", [True, False])
def test_chunking_partitions_the_weights(vds_lgt, raw_weights, split_multi):
    """One variant per chunk must give the same score as one chunk for all.

    Five weights rows, so `chunk_size=1` means five separate passes over the
    VDS, each with its own interval filter. If the chunk boundaries overlapped,
    a variant would be counted twice; if they left a gap, one would vanish.
    Either shows up here as a changed total.
    """
    n_single, single = totals(
        vds_lgt, raw_weights, chunk_size=None, split_multi=split_multi
    )
    n_per_variant, per_variant = totals(
        vds_lgt, raw_weights, chunk_size=1, split_multi=split_multi
    )

    assert n_single == 1
    assert n_per_variant == 5
    assert per_variant == single


def test_chunk_size_larger_than_the_table_is_a_single_chunk(
    vds_lgt, raw_weights
):
    """The default chunk_size (20000) dwarfs any test table; make sure the
    ceil() arithmetic doesn't produce a stray empty chunk."""
    n_chunks, _ = totals(vds_lgt, raw_weights, chunk_size=20000)
    assert n_chunks == 1


def test_aggregate_and_export_sums_across_chunks(
    vds_lgt, raw_weights, tmp_path, monkeypatch
):
    """End-to-end through the pandas aggregation and the CSV writer.

    `_aggregate_and_export` is the last step of `calculate_prs`, and the only
    part of the pipeline that needs GCS. Point `hfs.open` at a local file and
    the whole path is exercised.
    """
    out = tmp_path / "prs.csv"
    monkeypatch.setattr(
        _calculator.hfs, "open", lambda path, mode: open(path, mode)
    )

    config = PRSConfig(chunk_size=2)
    weights, n_chunks = _prepare_weights_for_chunking(
        weights_table=raw_weights, config=config, validate_table=True
    )
    partial_dfs = _process_chunks(
        full_weights_table=weights,
        n_chunks=n_chunks,
        vds=vds_lgt,
        config=config,
    )
    assert n_chunks == 3  # 5 variants, chunk_size 2 -> 2 + 2 + 1

    _aggregate_and_export(partial_dfs, str(out), config)

    written = pd.read_csv(out).set_index("person_id")["prs"].to_dict()
    # Same expectation as test_split_default_counts_alt_copies, now reassembled
    # from three separate passes over the VDS.
    assert written == {"S1": 0.0, "S2": 1.0, "S3": 3.0, "S4": 0.0}
