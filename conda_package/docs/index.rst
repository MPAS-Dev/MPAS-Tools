MPAS-Tools
==========

.. image:: images/so60to10.png
   :width: 500 px
   :align: center

MPAS-Tools includes a python package, compiled Fortran, C and C++ tools, and
scripts for supporting initialization, visualization and analysis of Model for
Prediction Across Scales (MPAS) components.  These tools are used by the
`COMPASS <https://github.com/MPAS-Dev/MPAS-Model/tree/master/testing_and_setup/compass>`_
(Configuring of MPAS Setups) framework within
`MPAS-Model <https://github.com/MPAS-Dev/MPAS-Model>`_ used to create
ocean and land-ice test cases,
the `MPAS-Analysis <https://github.com/MPAS-Dev/MPAS-Analysis>`_ package for
analyzing simulations, and in other MPAS-related workflows.

.. toctree::
   :caption: User's Guide
   :maxdepth: 2

   mesh_creation
   mesh_conversion
   interpolation

   cime

   config

   io

   logging

   transects

   vector

.. toctree::
   :caption: Visualization
   :maxdepth: 2

   visualization
   mpas_to_xdmf
   paraview_extractor

.. toctree::
   :caption: Landice Tools
   :maxdepth: 2

   landice/visualization

.. toctree::
   :caption: Ocean Tools
   :maxdepth: 2

   ocean/mesh_creation
   ocean/coastal_tools
   ocean/coastline_alteration
   ocean/moc
   ocean/depth
   ocean/streamfunction
   ocean/visualization

.. toctree::
   :caption: Sea-ice Tools
   :maxdepth: 2

   seaice/mask
   seaice/mesh
   seaice/partition
   seaice/regions
   seaice/regrid

.. toctree::
   :caption: Developer's Guide
   :maxdepth: 2

   making_changes
   testing_changes
   building_docs

   api

Indices and tables
==================

* :ref:`genindex`

.. toctree::
   :caption: Authors
   :maxdepth: 1

   authors
