Installation
============

You can install ``aoutools`` via ``pip`` using either the Python Package Index
(PyPI) or its GitHub repository. On the All of Us Researcher Workbench, you can
run the following commands directly in a Jupyter Notebook cell.

From PyPI
---------

Install the latest stable release from PyPI (recommended for production use):

.. code-block:: bash

   !pip install aoutools

From GitHub
-----------

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
