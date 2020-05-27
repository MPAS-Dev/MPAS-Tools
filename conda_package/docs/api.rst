#############
API reference
#############

This page provides an auto-generated summary of the MPAS mesh-tools API. For
more details and examples, refer to the relevant chapters in the main part of
the documentation.

MPAS mesh tools
===============

.. currentmodule:: mpas_tools.planar_hex

.. autosummary::
   :toctree: generated/

   make_planar_hex_mesh

.. currentmodule:: mpas_tools.translate

.. autosummary::
   :toctree: generated/

   translate

.. currentmodule:: mpas_tools.mesh.creation

.. autosummary::
   :toctree: generated/

   build_mesh.build_mesh
   coastal_tools.coastal_refined_mesh
   inject_bathymetry.inject_bathymetry
   inject_meshDensity.inject_meshDensity
   inject_preserve_floodplain.inject_preserve_floodplain
   jigsaw_driver.jigsaw_driver
   mesh_definition_tools.mergeCellWidthVsLat
   mesh_definition_tools.EC_CellWidthVsLat
   mesh_definition_tools.RRS_CellWidthVsLat
   mesh_definition_tools.AtlanticPacificGrid
   mpas_to_triangle.mpas_to_triangle
   open_msh.readmsh
   triangle_to_netcdf.triangle_to_netcdf
   jigsaw_to_netcdf.jigsaw_to_netcdf

.. currentmodule:: mpas_tools.mesh.conversion

.. autosummary::
   :toctree: generated/

   convert
   cull
   mask

.. currentmodule:: mpas_tools.merge_grids

.. autosummary::
   :toctree: generated/

   merge_grids

.. currentmodule:: mpas_tools.split_grids

.. autosummary::
   :toctree: generated/

   split_grids

.. currentmodule:: mpas_tools.io

.. autosummary::
   :toctree: generated/

   write_netcdf

.. currentmodule:: mpas_tools.scrip.from_mpas

.. autosummary::
   :toctree: generated/

   scrip_from_mpas

Ocean Tools
===========

.. currentmodule:: mpas_tools.ocean.coastline_alteration

.. autosummary::
   :toctree: generated/

   add_critical_land_blockages
   widen_transect_edge_masks

.. currentmodule:: mpas_tools.ocean.moc

.. autosummary::
   :toctree: generated/

   make_moc_basins_and_transects
   build_moc_basins

Visualization
=============

.. currentmodule:: mpas_tools.viz.paraview_extractor

.. autosummary::
   :toctree: generated/

   extract_vtk
