.. _dev_building_docs:

**************************
Building the Documentation
**************************

To make a local test build of the documentation, it is easiest to follow the
:ref:`dev_testing_changes` procedure for how to make a local build of the
``mpas_tools`` package.  The development environment includes the packages
needed to build the documentation. Simply run:

.. code-block:: bash

  cd conda_package/docs
  DOCS_VERSION=master make versioned-html

****************************
Previewing the Documentation
****************************

To preview the documentation locally, open the ``index.html`` file in the
``_build/html/master`` directory with your browser or try:

.. code-block:: bash

  cd _build/html
  python -m http.server 8000

Then, open http://0.0.0.0:8000/master/ in your browser.
