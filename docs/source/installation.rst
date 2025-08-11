Installation
============

You can install aoutools via pip using its GitHub repository.

On the All of Us Researcher Workbench, you can run this command directly in a
Jupyter Notebook cell:

.. code-block:: bash

    !pip install git+https://github.com/dokyoonkimlab/aoutools.git

aoutools has optional dependencies for the ``download_pgs`` function, which can
be used to download harmonized score files from `PGS Catalog`_. Note that
installing optional dependencies will disable the **dsub** job scheduler due to
a conflict with the pgscatalog.core package. The optional dependencies can be
installed via:

.. code-block:: bash

    !pip install "git+https://github.com/dokyoonkimlab/aoutools.git@main#egg=aoutools[pgs]"


.. _PGS Catalog: https://www.pgscatalog.org/
