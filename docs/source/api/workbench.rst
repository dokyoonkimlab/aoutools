:tocdepth: 3

aoutools (workbench)
====================

Helpers for setting up a Hail session on the *All of Us* Researcher Workbench.
Neither is required to use the rest of the library -- they exist so that the
Workbench-specific boilerplate (requester-pays billing, the reference genome,
and the location of the VariantDataset) lives in one place rather than being
copied into every notebook.

.. autofunction:: aoutools.init_hail

.. autofunction:: aoutools.get_vds_path

.. autodata:: aoutools.DEFAULT_VDS_PATH
   :annotation:
