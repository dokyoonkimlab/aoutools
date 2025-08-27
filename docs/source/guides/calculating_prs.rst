How-To: Calculate a Polygenic Risk Score (PRS)
==============================================

This guide shows you how to calculate one or more Polygenic Risk Scores (PRS)
using the `aoutools.prs` submodule.


.. note::

   To get started, you’ll need a Dataproc cluster set up for Hail Genomic
   Analysis. A good starting point is a master node with 8 CPUs, 52 GB RAM, and
   300 GB storage, plus 50 preemptible workers with 4 CPUs, 15 GB RAM, and 300
   GB storage each. This setup costs around $6 per hour and is enough to handle
   PRS weights from about 1 million variants. If you plan to run batch
   calculations with several large-scale weight files, you need a more powerful
   setup.


Step 1: Reading PRS Weights Files
---------------------------------
The `read_prs_weights` function is a flexible tool for importing weights files
into a validated Hail Table. It uses a `column_map` dictionary to handle
different file structures.

**Example 1: File with a header**

.. code-block:: python

    import hail as hl
    import logging
    from aoutools.prs import read_prs_weights

    # Set logging level to INFO to view logs and check function correctness.
    logging.basicConfig(level=logging.INFO)

    # Define a map from your file's column names to the required names
    column_map = {
        'chr': 'CHR',
        'pos': 'BP',
        'effect_allele': 'A1_EFFECT',
        'noneffect_allele': 'A2_NONEFFECT',
        'weight': 'BETA'
    }

    weights_ht = read_prs_weights(
        file_path='gs://my-bucket/data/my_weights.csv',
        header=True,
        column_map=column_map
    )

**Example 2: Header-less file (like PRS-CS output)**

The `column_map` uses 1-based integer indices instead of names. For convenience,
you can use the `read_prscs` wrapper if your weight file is generated from
PRScs.

.. code-block:: python

    from aoutools.prs import read_prscs

    # Using the main function
    column_map_no_header = {
        'chr': 1,
        'pos': 3,
        'effect_allele': 4,
        'noneffect_allele': 5,
        'weight': 6
    }

    weights_ht = read_prs_weights(
        file_path='local_weights_file.csv', # The function will auto-stage this to GCS
        header=False,
        column_map=column_map_no_header,
        keep_other_cols=True
    )


    # Using the convenient wrapper for PRS-CS files
    prscs_ht = read_prscs(
        file_path='local_prscs_output.csv'
    )


Step 2: Calculating a Single PRS
--------------------------------
Once you have a weights table and the All of Us VDS loaded, you can calculate a
PRS.

.. code-block:: python

    import os
    from aoutools.prs import calculate_prs, PRSConfig

    vds = hl.vds.read_vds(os.getenv('WGS_VDS_PATH'))

    # Assume 'weights_ht' is a Hail Table from Step 1
    prs_table = calculate_prs(
        weights_table=weights_ht,
        vds=vds,
        output_path='gs://my-bucket/results/single_prs.csv',
        config=PRSConfig()
    )

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

    from aoutools.prs import calculate_prs_batch

    # Create a dictionary mapping score names to their weights tables
    weights_map = {
        'CAD_prs': cad_weights_ht,
        'Asthma_prs': asthma_weights_ht,
    }

    # Calculate all scores in a single pass
    batch_prs_table = calculate_prs_batch(
        weights_map=weights_map,
        vds=vds,
        output_path='gs://my-bucket/results/batch_prs.csv',
        config=PRSConfig()
    )
