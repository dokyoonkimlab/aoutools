How-To: Using the ``download_pgs`` Function
===========================================

The ``download_pgs`` function is used to download harmonized score files from
the `PGS Catalog`_.

How It Works
------------

The ``download_pgs`` function is designed to resolve a specific dependency
conflict the can arise in the All of Us Researcher Workbench.

The core issue is a version mismatch between two important libraries:
``dsub``, which is used for job scheduling, and ``pgscatalog.core``, which is
used for downloading PGS data. They both rely on a shared dependency called
``tenacity``, but they require different, conflicting versions. Installing them
together would cause errors.

To solve this, ``download_pgs`` works as follows:

1. **Isolated Environment**: The first time you call the function, it
automatically creates a separate Python virtual environment using
``venv``. This keeps the dependencies for ``pgscatalog.core`` completely
isolated from your main notebook environment.

2. **Automatic Installation**: Inside this new, isolated environment, it
automatically installs ``pgscatalog.core`` and its required version of
``tenacity``.

3. **Execution**: It then runs the necessary download commands from within this
isolated, managed environment.

This method allows you to use both the ``dsub`` job scheduler and the
``download_pgs`` function in the same project without any dependency conflicts.

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
    from aoutools.prs import download_pgs

    # When this is called, the virtual environment will be created
    # at the custom path specified above.
    download_pgs(outdir='your_output_directory', pgs='PGS000001')


.. _PGS Catalog: https://www.pgscatalog.org/
