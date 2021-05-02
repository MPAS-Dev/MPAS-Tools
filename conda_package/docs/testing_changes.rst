.. _dev_testing_changes:

*****************************
Testing Changes to mpas_tools
*****************************

Here, we describe the workflow for creating a development conda environment
that points to ``mpas_tools`` in a branch from a local clone of the repo.
This approach works both for calling functions from the package within a python
script or another python package and for calling the "entry points"
(command-line tools; see :ref:`dev_making_changes`).

Basic instructions on how to install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
are beyond the scope of this documentation. Make sure the conda-forge channel
is added and that channel priority is "strict", meaning packages will
definitely come from conda-forge if they are available there.

.. code-block:: bash

    conda config --add channels conda-forge
    conda config --set channel_priority strict

To make a conda environment and install the current `mpas_tools` in a way that
it will be used out of the repo directly (i.e. it will notice changes as you
make them in your branch), run:

.. code-block:: bash

    cd conda_package
    conda env create -y -n mpas_tools_dev --file dev-spec.txt
    conda activate mpas_tools_dev
    python -m pip install -e .

You should now find that ``mpas_tools`` can be imported in python codes and the
various scripts and entry points are available in the path.

If you have already created the ``mpas_tools_dev`` environment, it may be best
to remove it (see below) and create it again.  If you are in a rush, you can
use:

.. code-block:: bash

    conda env update -f ./dev_environment
    conda activate mpas_tools_dev
    python -m pip install -e .

to update the existing environment and make sure ``mpas_tools`` in the
environment points to your current branch.

There is no need to build a conda package, as previous instructions had
suggested.

Removing the test environment
*****************************

If you're done with testing, you can remove the test environment

.. code-block::

  conda deactivate
  conda remove --all -n mpas_tools_dev
