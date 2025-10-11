How-To: Use the ``calculate_pgs`` Workflow
==========================================

The ``calculate_pgs`` function provides a streamlined, all-in-one workflow for
calculating Polygenic Risk Scores (PRS) directly from the PGS Catalog. It
automates the entire process, from downloading scoring files to calculating the
final scores.

This function is ideal when you know the specific PGS Catalog ID(s) you want to
analyze and prefer a single command to handle all the intermediate steps.

How It Works
------------

The ``calculate_pgs`` function simplifies the PRS calculation process by
combining three key steps into a single call:

1. **Download**: It begins by calling the ``download_pgs`` function internally
    to fetch the specified scoring files from the PGS Catalog and saves them to
    a temporary directory.
2. **Read**: It automatically reads each downloaded scoring file into a Hail
    Table, correctly mapping the standard PGS Catalog column names (e.g.,
    `hm_chr`, `hm_pos`, `effect_allele`, `other_allele`, `effect_weight`). If a
    file is malformed or cannot be read, it will be skipped with a warning.
3. **Calculate**: Finally, it uses the efficient ``calculate_prs_batch``
    function to calculate all requested PRS in a single pass, which minimizes
    reads of the Hail VDS.

Basic Usage
-----------

The following example demonstrates how to download two scores from the PGS
Catalog and calculate them for all samples in your VDS.

.. code-block:: python

    import os
    import pandas as pd
    import hail as hl

    # Import the workflow function
    from aoutools.prs import calculate_pgs

    # Initiate Hail and get bucket path
    hl.default_reference(new_default_reference="GRCh38")
    bucket = os.getenv("WORKSPACE_BUCKET")

    # Load the VDS
    vds = hl.vds.read_vds(os.getenv("WGS_VDS_PATH"))

    # Define the PGS IDs to calculate
    pgs_ids_to_calculate = ("PGS000196", "PGS000771")

    # Run the end-to-end workflow
    result_path = calculate_pgs(
        vds=vds,
        output_path=f"{bucket}/pgs_catalog_scores.csv",
        pgs=pgs_ids_to_calculate,
        build="GRCh38"
    )

    # Check the result
    pd.read_csv(result_path).head()

Important Considerations
------------------------

* **Input is Limited to PGS IDs**: To prevent accidental downloads of a large
  number of files, this function only accepts PGS Catalog IDs (e.g.,
  "PGS000123"). It does not support querying by EFO traits or PGP publication
  IDs.
* **Customization**: For advanced calculation options, you can pass a
  ``PRSConfig`` object to the ``config`` argument, just as you would with
  ``calculate_prs`` or ``calculate_prs_batch``.
* **Underlying Functions**: This workflow is a wrapper around other functions in
  the library. For more fine-grained control over downloading, reading, or
  calculation, you can use ``download_pgs``, ``read_prs_weights``, and
  ``calculate_prs_batch`` separately.
