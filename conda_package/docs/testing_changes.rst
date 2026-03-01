.. _dev_testing_changes:

*****************************
Testing Changes to mpas_tools
*****************************

Here, we describe the workflow for creating a development conda environment
that points to ``mpas_tools`` in a branch from a local clone of the repo.
The preferred workflow is to build and install the package locally with
``rattler-build`` from within a pixi environment.

.. _dev_local_rattler_build:

Building and Installing Locally with rattler-build
***************************************************

This is the recommended method, and it is required if you are modifying
compiled C++ or Fortran tools.

Use the pixi environment in ``conda_package/pixi.toml`` to run the build
command:

.. code-block:: bash

    cd conda_package
    pixi install
    pixi shell
    rattler-build build -m ci/linux_64_python3.14.____cpython.yaml -r recipe/ --output-dir ../output

This writes package artifacts to ``output/`` in the repository root.

To install the locally built package into an existing conda environment,
install the generated artifact file directly (replace the filename pattern with
the one produced on your platform):

.. code-block:: bash

    conda install -n mpas_tools_dev -y ../output/linux-64/mpas_tools-*.conda

If you use ``mamba`` or ``micromamba``, the same install command works with
those tools as well.

Quick-and-Dirty Alternative: Pixi Editable Install
**************************************************

If you are only making Python-level changes and do not need to rebuild the
compiled C++/Fortran tools, you can use editable installation in pixi:

.. code-block:: bash

    cd conda_package
    pixi install
    pixi shell
    pixi run install-editable

Then run tools within the pixi shell (for example ``pytest``).

.. important::

    Editable installation updates Python code but does **not** rebuild compiled
    C++ and Fortran command-line tools.

    A useful hybrid workflow is to install the latest release conda package
    first (to get compiled tools), then install your branch in editable mode on
    top for Python development.

Legacy Method: Conda Editable Install
*************************************

This workflow is kept for compatibility but is no longer the preferred method.

Basic instructions on how to install
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ are beyond the
scope of this documentation. Make sure the conda-forge channel is added and
that channel priority is "strict":

.. code-block:: bash

    conda config --add channels conda-forge
    conda config --set channel_priority strict

Then create and activate a development environment from ``dev-spec.txt`` and
install in editable mode:

.. code-block:: bash

    cd conda_package
    conda env create -y -n mpas_tools_dev --file dev-spec.txt
    conda activate mpas_tools_dev
    python -m pip install --no-deps --no-build-isolation -e .

Removing the test environment
-----------------------------

If you're done with testing, you can remove the test environment

.. code-block::

  conda deactivate
  conda remove --all -n mpas_tools_dev
