> **Disclaimer:** This project is **not affiliated with, endorsed by, or
> sponsored by** the *All of Us Research Program*. The software is provided **as
> is**, without warranty. It is in an **early development stage**, and
> functions, APIs, and signatures may change frequently without notice.

# aoutools: Tools for All of Us Researcher Workbench

aoutools is a Python library designed to simplify common analysis tasks on the
All of Us Researcher Workbench.

## Overview

aoutools aims to provides a suite of high-level functions designed to make All
of Us data analyses more accessible.

The initial release focuses on the `aoutools.prs` submodule, which offers
convenient functions for:

1.  **Reading PRS Weights Files:** A flexible reader that can handle various
    file formats, with and without headers.
2.  **Calculating PRS:** A cost-efficient strategy for calculating PRS directly
    on the All of Us VDS, with support for batch mode to calculate multiple
    scores at once.

## Installation

You can install aoutools via pip on All of Us Workbench:

```bash
# Use !pip if running in a Jupyter notebook
pip install git+https://github.com/dokyoonkimlab/aoutools.git
```

Please check the online documentation at: [aoutools
Documentation](https://aoutools.readthedocs.io)

## Usage

All functions are available under the `aoutools.prs` submodule.

### 1. Reading PRS Weights Files

The `read_prs_weights` function is a flexible tool for importing weights files
into a validated Hail Table. It uses a `column_map` dictionary to handle
different file structures.

#### Example 1: Reading a file with a header

```python
import hail as hl
from aoutools.prs import read_prs_weights

# Define a map from your file's column names to the required names
column_map = {
    'chr': 'CHR',
    'pos': 'BP',
    'effect_allele': 'A1_EFFECT',
    'noneffect_allele': 'A2_NONEFFECT',
    'weight': 'BETA'
}

weights_ht = read_prs_weights(
    file_path='gs://my-bucket/data/my_weights.tsv',
    header=True,
    column_map=column_map,
    keep_other_cols=True  # Optional: keeps other columns from your file
)

weights_ht.show()
```

#### Example 2: Reading a header-less file (like PRS-CS output)

For files without a header, the `column_map` uses 1-based integer indices
instead of names.

```python
from aoutools.prs import read_prs_weights, read_prscs

# Using the main function
column_map_no_header = {
    'chr': 1,
    'pos': 3,
    'effect_allele': 4,
    'noneffect_allele': 5,
    'weight': 6
}

weights_ht = read_prs_weights(
    file_path='local_weights_file.txt', # The function will auto-stage this to GCS
    header=False,
    column_map=column_map_no_header,
    keep_other_cols=True
)

# Using the convenient wrapper for PRS-CS files
prscs_ht = read_prscs(
    file_path='local_prscs_output.txt',
    keep_other_cols=True
)
```

### 2. Calculating a Polygenic Risk Score

Once you have a weights table and have loaded the All of Us VDS, you can
calculate a PRS.

```python
import os
from aoutools.prs import calculate_prs, PRSConfig

# Load the All of Us VDS
vds_path = os.getenv('WGS_VDS_PATH')
vds = hl.vds.read_vds(vds_path)

# Assume 'weights_ht' is a Hail Table from read_prs_weights
config = PRSConfig()

prs_table = calculate_prs(
    weights_table=weights_ht,
    vds=vds,
    output_path='gs://my-bucket/results/single_prs.tsv',
    config=config
)

# You can then read the result file
if output_file_path:
    import pandas as pd
    df = pd.read_csv(output_file_path, sep='\t')
    print(df.head())
```

#### Handling Odds Ratios (OR)

If your weights file contains Odds Ratios instead of BETA values, the function
can automatically log-transform them.

```python
# Assume 'or_weights_ht' has a column named 'OR'
config_or = PRSConfig(
    weight_col_name='OR',      # Specify the weight column
    log_transform_weight=True  # Enable log transformation
)

calculate_prs(
    weights_table=or_weights_ht,
    vds=vds,
    output_path='gs://my-bucket/results/or_prs.tsv',
    config=config_or
)
```

### 3. Batch PRS Calculation

To calculate multiple scores efficiently, use the `calculate_prs_batch`
function. This is highly recommended as it reads the VDS only once.

```python
from aoutools.prs import calculate_prs_batch

# Create a dictionary mapping score names to their weights tables
weights_map = {
    'CAD_prs': cad_weights_ht,
    'Asthma_prs': asthma_weights_ht,
    'T2D_prs': t2d_weights_ht
}

# Calculate all scores in a single pass
batch_prs_table = calculate_prs_batch(
    weights_map=weights_map,
    vds=vds
)

# You can then use this path to read the results.
if output_file_path:
    import pandas as pd
    df = pd.read_csv(output_file_path, sep='\t')
    print(df.head())
```
