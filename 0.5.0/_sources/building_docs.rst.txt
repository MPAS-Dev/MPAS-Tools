.. _dev_building_docs:

**************************
Building the Documentation
**************************

To make a local test build of the documentation, it is easiest to follow the
:ref:`dev_testing_changes` procedure for how to make a local build of the
``mpas_tools`` package.  Then, you need to set up a conda environment with the
test build and some other required packages:

code-block::

  $ conda create -y -n test_mpas_tools_docs --use-local mpas_tools sphinx mock \
       sphinx_rtd_theme
  $ conda activate test_mpas_tools_docs

Then, to build the documentation, run:

code-block::

  $ export DOCS_VERSION="test"
  $ cd conda_package/docs
  $ make html

Then, you can view the documentation by opening ``_build/html/index.html``.
