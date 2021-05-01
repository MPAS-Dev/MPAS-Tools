.. _dev_building_docs:

**************************
Building the Documentation
**************************

To make a local test build of the documentation, it is easiest to follow the
:ref:`dev_testing_changes` procedure for how to make a local build of the
``mpas_tools`` package.  The development environment includes the packages
needed to build the documentation. Simply run:

code-block::

  export DOCS_VERSION="test"
  cd conda_package/docs
  make html

Then, you can view the documentation by opening ``_build/html/index.html``.
