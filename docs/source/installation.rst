Installation
============

You can install ``aoutools`` via pip using its GitHub repository. On the All of
Us Researcher Workbench, you can run this command directly in a Jupyter Notebook
cell:

.. code-block:: bash

    !pip install git+https://github.com/dokyoonkimlab/aoutools.git

Using the ``download_pgs`` Function
===================================

The ``download_pgs`` function is used to download harmonized score files from
the `PGS Catalog`_.

How It Works
------------

The design of the ``download_pgs`` function is centered on resolving a specific
dependency conflict to ensure compatibility within the All of Us Researcher
Workbench.

The core issue is that the Google Cloud **dsub** library and the
**``pgscatalog.core``** library require conflicting versions of a shared
dependency, ``tenacity``. Installing both in the same environment would cause
errors.

To solve this, ``download_pgs`` works as follows:

1. **Isolated Environment**: Upon its first use, the function automatically
    creates a separate Python virtual environment using ``venv``. This keeps the
    dependencies for ``pgscatalog.core`` completely separate from your main
    notebook environment.
2. **Automatic Installation**: Inside this isolated environment, it installs
    ``pgscatalog.core`` and its specific version of ``tenacity``.
3. **Execution**: It then runs the required download commands from within this
    managed environment.

This approach allows you to use both the ``dsub`` job scheduler and the
``download_pgs`` function in the same project without any conflicts.

Customizing the Virtual Environment Path
----------------------------------------

By default, the virtual environment is created in your home directory at
``~/.aoutools/pgscatalog_env``.

You can override this location by setting the ``AOUTOOLS_PGS_ENV_DIR``
environment variable. This must be set *before* the ``download_pgs`` function is
called for the first time. To set this in a Jupyter Notebook, you can use the
``os`` module.

For example, to create the environment in a different directory, add the
following to a notebook cell:

.. code-block:: python

    import os

    # Set the environment variable to a custom path before importing aoutools
    os.environ['AOUTOOLS_PGS_ENV_DIR'] = '/path/to/your/custom/env/dir'

    # Now you can import and use the function
    from aoutools.pgs import download_pgs

    # When this is called, the virtual environment will be created
    # at the custom path specified above.
    download_pgs(outdir='your_output_directory', pgs='PGS000001')

.. _PGS Catalog: https://www.pgscatalog.org/

