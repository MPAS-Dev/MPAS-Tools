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
    rattler-build build -m ci/linux_64_python3.14.____cpython.yaml -r recipe/ --output-dir output

This writes package artifacts to ``output/`` under ``conda_package``.

To install the locally built package into the pixi environment, add the local
build output as a channel and then add ``mpas_tools`` from that channel:

.. code-block:: bash

    cd conda_package
    pixi workspace channel add "file://$PWD/output"
    pixi add --platform linux-64 "mpas_tools [channel='file://$PWD/output']"

.. important::

    pixi, like other conda-family package managers, identifies a package by
    its name, version and build string.  If you rebuild a local package without
    changing that identity, pixi may continue using an older cached artifact
    even if the file in ``output/`` has changed.

    If you rebuild ``mpas_tools`` locally and need pixi to pick up the new
    package contents reliably, bump the conda recipe build number in
    ``conda_package/recipe/recipe.yaml`` before rebuilding.  For Python-only
    development, ``pixi run install-editable`` is often more convenient.

.. warning::

    The commands above modify ``pixi.toml`` (and possibly ``pixi.lock``).
    These are local development changes only and should **not** be committed.
    Before opening a PR, reset those files to the repository state.

On macOS, use ``--platform osx-64`` instead.  If your workspace includes both
``linux-64`` and ``osx-64`` in ``pixi.toml``, pixi may try to solve both
platforms during dependency updates.  In that case, build local artifacts for
both platforms (with corresponding CI recipe files) or add the local package
for only the platform you built.

If you want to return to using only published channels afterward, you can
remove the local channel from ``pixi.toml``.

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
    top for Python development.  This also avoids the need to bump the conda
    build number for every local Python-only rebuild.

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
