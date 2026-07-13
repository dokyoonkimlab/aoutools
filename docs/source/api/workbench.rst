:tocdepth: 3

aoutools (workbench)
====================

Helpers for setting up a Hail session on the *All of Us* Researcher Workbench.
Neither is required to use the rest of the library -- they exist so that the
Workbench-specific boilerplate (requester-pays billing, the reference genome,
and the location of the VariantDataset) lives in one place rather than being
copied into every notebook.

.. note::

   The current *All of Us* platform (Verily) no longer exports the
   ``WGS_VDS_PATH`` and ``WORKSPACE_BUCKET`` environment variables that the
   older Workbench guaranteed. Code that reads them directly fails on a fresh
   workspace. These helpers resolve both without them.

   Exporting a variable from a Jupyter **terminal** does not reach the
   notebook: the terminal and the kernel are sibling children of the Jupyter
   server, which captured its environment at startup, so the kernel never sees
   it -- even after a restart. Set it from inside a notebook cell instead, or
   put it in ``~/.ipython/profile_default/startup/00-aou-env.py``.

.. autofunction:: aoutools.init_hail

.. autofunction:: aoutools.get_vds_path

.. autofunction:: aoutools.get_workspace_bucket

.. autodata:: aoutools.DEFAULT_VDS_PATH
   :annotation:
