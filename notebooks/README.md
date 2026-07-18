# Validation notebooks

These notebooks are the **human-run validation tier**. They are *not* run by CI —
someone runs them by hand before a release to check the library against the real
*All of Us* data, on holes the automated `tests/` suites can't reach. A few also
serve as **evidence records**: they document why a known issue is closed, so
re-running one re-confirms the finding.

Most need a real VDS and so run only on the Workbench (a **Hail Genomic Analysis**
environment). One runs offline. See the last column.

| Notebook | What it checks | Where it runs |
|---|---|---|
| `validate_scoring_on_aou.ipynb` | Finds a real variant of each shape (biallelic, multi-allelic, non-minimal, no-call) and checks scores against an independent oracle (`hl.vds.to_dense_mt`). The real-data counterpart of `tests/integration/`. | Workbench |
| `validate_public_api_on_aou.ipynb` | The user-facing functions and the `gs://` output path — `calculate_prs`, `calculate_prs_batch`, `calculate_pgs`, the PGS Catalog download, local→bucket staging — none of which any offline test reaches. | Workbench |
| `validate_synthetic_control_on_aou.ipynb` | Builds a synthetic dataset with a known answer and drives the public API against it, including a per-variant diagnostic toolkit for any discrepancy. | Offline (`pixi run -e integration`), except the `gs://` cells → Workbench |
| `verify_split_multi_downcoding.ipynb` | Confirms that splitting a multi-allelic site relabels other ALTs as reference — the fact that makes two dosages necessary. | Workbench |
| `measure_minrep_locus_shift.ipynb` | Measures the locus-shift rate in the real VDS (it is zero); re-run if the scoring tripwire ever fires. | Workbench |

## Running one

On the Workbench, start a **Hail Genomic Analysis** environment, then in the first
cell install the library and bootstrap Hail:

```python
# !pip install -q git+https://github.com/dokyoonkimlab/aoutools.git@dev
from aoutools import init_hail
init_hail()
```

`init_hail()` wires up requester-pays billing and the GRCh38 reference; the
notebooks resolve the VDS path themselves via `get_vds_path()`.

The offline cells of `validate_synthetic_control_on_aou.ipynb` need no Workbench —
run them locally under `pixi run -e integration`.

## Benchmark (not a validation tier)

`benchmark_scoring_speed.ipynb` times `calculate_prs` on a real VDS window. It
asserts nothing about correctness -- run it once per install to compare scoring
speed across versions (e.g. before merging a performance change). Workbench only.

## Why each notebook exists

For the design reasoning behind these — which correctness finding each one pins,
and why the automated tiers can't cover it — see the **"third tier"** section of
`AGENTS.md`.
