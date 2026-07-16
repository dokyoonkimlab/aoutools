How-To: Calculate a Polygenic Risk Score (PRS)
==============================================

This guide shows you how to calculate one or more Polygenic Risk Scores (PRS)
using the **aoutools.prs** submodule.


.. note::

   To get started, you’ll need a cluster configured for Hail Genomic Analysis.
   A reasonable starting point is an ``n2-standard-8`` (or larger) main node
   together with 30 secondary workers running as Spot VMs, each keeping the
   default 500 GB disk. This handles PRS weights from about 1 million variants
   across all available All of Us WGS samples, with no variant or sample
   filtering. Spot VMs keep the cost down but can be reclaimed by the platform
   at any time; for a run that finishes in a couple of minutes, the chance of an
   interruption is low. Exact pricing depends on your region and current Spot
   rates, but this configuration comes to roughly a few US dollars per hour. If
   you plan to run batch calculations across several large-scale weight files,
   choose a more powerful setup.

Given the specified workspace configuration, the tutorial requires only 1–2
minutes for both single and batch PRS calculations, not including the
file-reading time, which is under a minute.

Setup
-----

.. code-block:: python

    import logging
    from importlib import resources
    import pandas as pd
    import hail as hl

    # Workbench session helpers
    from aoutools import init_hail, get_vds_path, get_workspace_bucket

    # Import functions to calculate PRS
    from aoutools.prs import (
        read_prs_weights,
        calculate_prs,
        calculate_prs_batch,
        PRSConfig
    )

    # Set logging level to INFO to view logs and check function correctness.
    logging.basicConfig(level=logging.INFO)

    # Initialize Hail for the Workbench. This wires up requester-pays billing
    # and sets the GRCh38 reference for you.
    init_hail()

    # Workspace bucket for output files
    bucket = get_workspace_bucket()

    # Get paths for example data
    data_dir_path = str(resources.files("aoutools.data"))
    prs_weights_header = f"{data_dir_path}/prs_weights_header.csv"
    prs_weights_noheader = f"{data_dir_path}/prs_weights_noheader.tsv"

    # Load the All of Us WGS VariantDataset
    vds = hl.vds.read_vds(get_vds_path())


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

**Example 2: Header-less file**

For a header-less file, the `column_map` uses 1-based integer indices instead of
names.

.. code-block:: python

    column_map_noheader = {
        'chr': 1,
        'pos': 3,
        'effect_allele': 4,
        'noneffect_allele': 5,
        'weight': 6
    }

    weights_ht_noheader = read_prs_weights(
        file_path=prs_weights_noheader,
        header=False,
        column_map=column_map_noheader,
        delimiter='\t'
    )


Step 2: Calculating a Single PRS
--------------------------------
Once you have a weights table and the All of Us VDS loaded, you can calculate a
PRS. To see more detailed log information, set
``PRSConfig(detailed_timings=True)`` and pass it to the ``config`` argument.

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
To calculate multiple scores efficiently, use ``calculate_prs_batch``. This is
highly recommended as it reads the VDS only once.

.. code-block:: python

    # Create a dictionary mapping score names to their weights tables
    weights_tables_map = {
        'prs1': weights_ht_header,
        'prs2': weights_ht_noheader,
    }

    # Calculate all scores in a single pass
    prs_batch = calculate_prs_batch(
        weights_tables_map=weights_tables_map,
        vds=vds,
        output_path=f"{bucket}/batch_prs.csv"
    )

    # Check the result
    pd.read_csv(prs_batch).head()
