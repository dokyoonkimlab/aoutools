How-To: Calculate a Polygenic Risk Score (PRS)
==============================================

This guide shows you how to calculate one or more Polygenic Risk Scores (PRS)
using the **aoutools.prs** submodule.


.. note::

   To get started, you’ll need a Dataproc cluster set up for Hail Genomic
   Analysis. A good starting point is a master node with 8 CPUs, 52 GB RAM, and
   300 GB storage, plus 50 preemptible workers with 4 CPUs, 15 GB RAM, and 300
   GB storage each. This setup costs around $6 per hour and is enough to handle
   PRS weights from about 1 million variants. If you plan to run batch
   calculations with several large-scale weight files, you need a more powerful
   setup.


Setup
-----

.. code-block:: python

    import os
    import logging
    from importlib import resources
    import pandas as pd
    import hail as hl

    # Import functions to calculate PRS
    from aoutools.prs import (
        read_prs_weights,
        read_prscs,
        calculate_prs,
        calculate_prs_batch,
        PRSConfig
    )

    # Set logging level to INFO to view logs and check function correctness.
    logging.basicConfig(level=logging.INFO)

    # Bucket path
    bucket = os.getenv("WORKSPACE_BUCKET")

    # Get paths for example data
    data_dir_path = str(resources.files("aoutools.data"))
    prs_weights_header = f"{data_dir_path}/prs_weights_header.csv"
    prs_weights_noheader = f"{data_dir_path}/prs_weights_noheader.tsv"

    # Initiate Hail
    hl.default_reference(new_default_reference="GRCh38")

    # VDS file
    vds = hl.vds.read_vds(os.getenv("WGS_VDS_PATH"))


Step 1: Reading PRS Weights Files
---------------------------------
The ``read_prs_weights`` function is a flexible tool for importing weights files
into a validated Hail Table. It uses a `column_map` dictionary to handle
different file structures.

**Example 1: File with a header**

.. code-block:: python

    # Define a map from your file's column names to the required names
    column_map_header = {
        'chr': 'CHR',
        'pos': 'POS',
        'effect_allele': 'A1',
        'noneffect_allele': 'A2',
        'weight': 'WEIGHT'
    }

    weights_ht_header = read_prs_weights(
        file_path=prs_weights_header,
        header=True,
        column_map=column_map_header
    )

.. note::

   The files `prs_weights_header` and `prs_weights_noheader` are local, but
   Hail cannot access them directly from a local Jupyter environment. Therefore,
   the ``read_prs_weights`` function automatically stages an input file to a
   temporary Google Cloud Storage (GCS) location at
   gs://your-workspace-bucket/data/temp_prs_data so that Hail can access them.
   However, it is recommended for users to upload the input files to a GCS
   bucket and provide a path that starts with 'gs://'.

**Example 2: Header-less file (like PRScs output)**

The `column_map` uses 1-based integer indices instead of names. For convenience,
you can use the ``read_prscs`` wrapper if your weight file is generated from
PRScs.

.. code-block:: python

    # Using the main function
    column_map_noheader = {
        'chr': 1,
        'pos': 3,
        'effect_allele': 4,
        'noneffect_allele': 5,
        'weight': 6
    }

    weights_ht_noheader = read_prs_weights(
        file_path=prs_weights_noheader
        header=False,
        column_map=column_map_noheader
    )


    # Using the convenient wrapper for PRS-CS files
    prscs_ht = read_prscs(
        file_path=prs_weights_noheader
    )


Step 2: Calculating a Single PRS
--------------------------------
Once you have a weights table and the All of Us VDS loaded, you can calculate a
PRS.

.. code-block:: python

    # Assume 'weights_ht_header' is a Hail Table from Step 1
    prs_single = calculate_prs(
        weights_table=weights_ht_header,
        vds=vds,
        output_path=f"{bucket}/single_prs.csv"
    )

    # Check the result
    pd.read_csv(prs_single).head()


**Advanced: Handling Odds Ratios (OR)**

If your weights file uses Odds Ratios, the function can log-transform them into
BETA values.

.. code-block:: python

    config_or = PRSConfig(
        weight_col_name='OR',
        log_transform_weight=True
    )


Tip: Batch PRS Calculation
------------------------------------
To calculate multiple scores efficiently, use `calculate_prs_batch`. This is
highly recommended as it reads the VDS only once.

.. code-block:: python

    # Create a dictionary mapping score names to their weights tables
    weights_map = {
        'prs1': weights_ht_header,
        'prs2': weights_ht_noheader,
    }

    # Calculate all scores in a single pass
    prs_batch = calculate_prs_batch(
        weights_map=weights_map,
        vds=vds,
        output_path=f"{bucket}/batch_prs.csv"
    )

    # Check the result
    pd.read_csv(prs_batch).head()
