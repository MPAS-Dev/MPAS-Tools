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

   build_mesh.build_spherical_mesh
   build_mesh.build_planar_mesh
   jigsaw_driver.jigsaw_driver
   jigsaw_to_netcdf.jigsaw_to_netcdf
   mesh_definition_tools.mergeCellWidthVsLat
   mesh_definition_tools.EC_CellWidthVsLat
   mesh_definition_tools.RRS_CellWidthVsLat
   mesh_definition_tools.AtlanticPacificGrid
   mpas_to_triangle.mpas_to_triangle
   open_msh.readmsh
   signed_distance.signed_distance_from_geojson
   signed_distance.mask_from_geojson
   signed_distance.distance_from_geojson
   triangle_to_netcdf.triangle_to_netcdf

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

.. currentmodule:: mpas_tools.ocean.coastal_tools

.. autosummary::
   :toctree: generated/

   coastal_refined_mesh
   create_background_mesh
   extract_coastlines
   distance_to_coast
   compute_cell_width
   save_matfile
   CPP_projection
   smooth_coastline
   get_data_inside_box
   get_indices_inside_quad
   get_convex_hull_coordinates
   plot_coarse_coast
   plot_region_box

.. currentmodule:: mpas_tools.ocean.inject_bathymetry

.. autosummary::
   :toctree: generated/

   inject_bathymetry

.. currentmodule:: mpas_tools.ocean.inject_meshDensity

.. autosummary::
   :toctree: generated/

   inject_meshDensity

.. currentmodule:: mpas_tools.ocean.inject_preserve_floodplain

.. autosummary::
   :toctree: generated/

   inject_preserve_floodplain


Visualization
=============

.. currentmodule:: mpas_tools.viz.paraview_extractor

.. autosummary::
   :toctree: generated/

   extract_vtk
