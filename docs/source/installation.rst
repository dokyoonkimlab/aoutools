Installation
============

Requirements
------------

``aoutools`` targets the **All of Us Researcher Workbench 2.0**. Version 0.2.0
is the first release for it: it adds :func:`~aoutools.init_hail` and
:func:`~aoutools.get_vds_path`, which wire up the new VDS location and the Hail
setup that platform expects. Releases 0.1.x targeted Workbench 1.0, which was
decommissioned on June 30, 2026.

The PRS functions read the *All of Us* VDS through Hail, so they need a **Hail
Genomic Analysis** environment rather than a general analysis environment. Hail
itself comes with the Workbench and ``aoutools`` deliberately does not install
it: Hail pins ``numpy`` and ``pandas`` tightly, and installing a second copy
could disturb the versions your environment already provides.

Python 3.11 or newer is required. The Workbench genomics runtime is Python 3.11.

Creating the cloud environment
------------------------------

Hail distributes its work across a Spark cluster, so the cloud environment must
be configured before ``aoutools`` is installed. In the Workbench's cloud
environment panel, select the following:

Environment type
   **JupyterLab Spark cluster**. The PRS functions cannot run in a standard
   (non-Spark) environment.

Software to install
   **Hail (Spark 3.5.3, hail 0.2.135)**. These are the exact versions
   ``aoutools`` is tested against — the continuous integration environment pins
   the same ones — so selecting a different combination is untested.

Master node
   **n2-standard-8** or larger.

Secondary workers
   **Typically 10 to 50**, and more for a large job. Adding workers shortens the
   run, but the speed-up is not proportional to the number added, so choose a
   count based on the size of the score you are calculating and how soon you
   need the result rather than assuming more is always better.

Once the environment is running, install ``aoutools`` as below.

Installing the package
----------------------

You can install ``aoutools`` via ``pip`` using either the Python Package Index
(PyPI) or its GitHub repository. On the All of Us Researcher Workbench, you can
run the following commands directly in a Jupyter Notebook cell.

From PyPI
~~~~~~~~~

Install the latest stable release from PyPI (recommended for production use):

.. code-block:: bash

   !pip install aoutools

From GitHub
~~~~~~~~~~~

Install the latest version from the main branch on GitHub. This may include new
features or bug fixes not yet released on PyPI.

.. code-block:: bash

   !pip install git+https://github.com/dokyoonkimlab/aoutools.git


Troubleshooting
---------------

If the package installs successfully but you encounter a ``ModuleNotFoundError``
when trying to import it, please restart the Jupyter Notebook kernel. This can
happen when a new package is installed in the current environment, and
restarting the kernel ensures that the new package is properly loaded and
recognized.
