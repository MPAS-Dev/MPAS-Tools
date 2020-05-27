.. _dev_testing_changes:

*****************************
Testing Changes to mpas_tools
*****************************

There are a few different ways to test the ``mpas_tools`` package.  Typically,
the quickest turn-around between making changes and seeing their results are
going to be seen if you can test the code straight out of the git repo.  This
approach works for calling functions from the package within a python script
but doesn't give you easy access to the "entry points"(see
:ref:`dev_making_changes`). To more fully test the package, you will need to
build the package locally, install it into a new conda environment, and test
your code within that environment.

Testing from the git repo
=========================

If you are testing a simple python script that accesses ``mpas_tools``, you can
make a symlink in the same directory as your python script to ``mpas_tools``
within ``conda_package``.  Python should search the local path before looking
elsewhere so this should work even if a previous version ``mpas_tools`` is
already installed in the conda environment you are using.

Testing the conda package
=========================

Updating the Version
********************

As part of your testing, you should update the version of ``mpas_tools``.  This
should be done both in ``conda_package/mpas_tools/__init__.py``:

.. code-block:: python

  __version_info__ = (0, 0, 11)
  __version__ = '.'.join(str(vi) for vi in __version_info__)

Increment ``__version_info__`` (major, minor or micro version, depending on
what makes sense).

The version in the conda recipe (``conda_package/recipe/meta.yaml``) needs to
match:

.. code-block::

  {% set name = "mpas_tools" %}
  {% set version = "0.0.11" %}

It is also a good idea to add the new version to the :ref:`versions`.  The new
links won't be valid until a new release is made and Travis CI has generated
the associated documentation.  Eventually, it should be possible to do this
automatically but that has not yet been implemented.

Building the package
********************

To build the conda package, you will need to install conda-build into your base
conda environment.  (Basic instructions on how to install Miniconda or Anaconda
are beyond the scope of this documentation.)

.. code-block::

  $ conda config --add channels conda-forge
  $ conda config --set channel_priority strict
  $ conda install -n base conda-build

To build the package, make sure you are in the base of the repo and run:

.. code-block::

  $ rm -rf ~/miniconda3/conda-bld
  $ conda build conda_package/recipe

The first is to make sure you don't have existing packages already built that
would get used in your building and testing instead of the versions from
``conda-forge``.  If your conda setup is installed somewhere other than
``~/miniconda3``, use the appropriate path.

Installing the package
**********************

To make a new test environment to try out scripts, other python packages or
other workflows that use the tools, run:

.. code-block::

  $ conda create -n test_mpas_tools --use-local python=3.8 mpas_tools

You can name the environment whatever if useful to you.  Activate the
environment with:

.. code-block::

  $ conda activate test_mpas_tools

You should now find that ``mpas_tools`` can be imported in python codes and the
various scripts and entry points are available in the path.

Removing the test environment
*****************************

If you're done with testing, you can remove the test environment

.. code-block::

  $ conda deactivate
  $ conda remove --all -n test_mpas_tools

