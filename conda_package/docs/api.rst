#############
API reference
#############

This page provides an auto-generated summary of the MPAS mesh-tools API. For
more details and examples, refer to the relevant chapters in the main part of
the documentation.

MPAS mesh tools
===============

Mesh creation
-------------

.. currentmodule:: mpas_tools.planar_hex

.. autosummary::
   :toctree: generated/

   make_planar_hex_mesh

.. currentmodule:: mpas_tools.mesh.creation

.. autosummary::
   :toctree: generated/

   build_mesh
   build_mesh.build_spherical_mesh
   build_mesh.build_planar_mesh
   jigsaw_driver.jigsaw_driver
   jigsaw_driver.build_jigsaw
   jigsaw_to_netcdf.jigsaw_to_netcdf
   mesh_definition_tools
   mesh_definition_tools.mergeCellWidthVsLat
   mesh_definition_tools.EC_CellWidthVsLat
   mesh_definition_tools.RRS_CellWidthVsLat
   mesh_definition_tools.AtlanticPacificGrid
   mpas_to_triangle.mpas_to_triangle
   signed_distance
   signed_distance.signed_distance_from_geojson
   signed_distance.mask_from_geojson
   signed_distance.distance_from_geojson
   triangle_to_netcdf.triangle_to_netcdf

Mesh conversion
---------------

.. currentmodule:: mpas_tools.mesh.conversion

.. autosummary::
   :toctree: generated/

   convert
   cull
   mask

.. currentmodule:: mpas_tools.mesh.cull

.. autosummary::
   :toctree: generated/

   write_map_culled_to_base
   map_culled_to_base
   write_culled_dataset
   cull_dataset

.. currentmodule:: mpas_tools.mesh.mask

.. autosummary::
   :toctree: generated/

   compute_mpas_region_masks
   compute_mpas_transect_masks
   compute_mpas_flood_fill_mask
   compute_lon_lat_region_masks

.. currentmodule:: mpas_tools.merge_grids

.. autosummary::
   :toctree: generated/

   merge_grids

.. currentmodule:: mpas_tools.split_grids

.. autosummary::
   :toctree: generated/

   split_grids

.. currentmodule:: mpas_tools.translate

.. autosummary::
   :toctree: generated/

   translate
   center
   center_on_mesh

.. currentmodule:: mpas_tools.scrip.from_mpas

.. autosummary::
   :toctree: generated/

   scrip_from_mpas

Config
------

.. currentmodule:: mpas_tools.config

.. autosummary::
   :toctree: generated/

   MpasConfigParser
   MpasConfigParser.add_user_config
   MpasConfigParser.add_from_file
   MpasConfigParser.add_from_package
   MpasConfigParser.get
   MpasConfigParser.getint
   MpasConfigParser.getfloat
   MpasConfigParser.getboolean
   MpasConfigParser.getlist
   MpasConfigParser.getexpression
   MpasConfigParser.has_section
   MpasConfigParser.has_option
   MpasConfigParser.set
   MpasConfigParser.write
   MpasConfigParser.copy
   MpasConfigParser.append
   MpasConfigParser.prepend
   MpasConfigParser.__getitem__

I/O
---

.. currentmodule:: mpas_tools.io

.. autosummary::
   :toctree: generated/

   write_netcdf

Parallelism
-----------

.. currentmodule:: mpas_tools.parallel

.. autosummary::
   :toctree: generated/

   create_pool

Interpolation
-------------

.. currentmodule:: mpas_tools.mesh.interpolation

.. autosummary::
   :toctree: generated/

   interp_bilin

CIME constants
--------------

.. currentmodule:: mpas_tools.cime

.. autosummary::
   :toctree: generated/

   constants

Landice Tools
=============

.. currentmodule:: mpas_tools.landice

.. autosummary::
   :toctree: generated/

    visualization
    visualization.plot_transect
    visualization.plot_map
    visualization.plot_grounding_lines

Ocean Tools
===========

.. currentmodule:: mpas_tools.ocean

.. autosummary::
   :toctree: generated/

   coastline_alteration
   coastline_alteration.add_critical_land_blockages
   coastline_alteration.widen_transect_edge_masks
   coastline_alteration.add_land_locked_cells_to_mask

   moc
   moc.make_moc_basins_and_transects
   moc.add_moc_southern_boundary_transects

   build_mesh
   build_mesh.build_spherical_mesh
   build_mesh.build_planar_mesh

   coastal_tools
   coastal_tools.coastal_refined_mesh
   coastal_tools.create_background_mesh
   coastal_tools.extract_coastlines
   coastal_tools.distance_to_coast
   coastal_tools.compute_cell_width
   coastal_tools.save_matfile
   coastal_tools.CPP_projection
   coastal_tools.smooth_coastline
   coastal_tools.get_data_inside_box
   coastal_tools.get_indices_inside_quad
   coastal_tools.get_convex_hull_coordinates
   coastal_tools.plot_coarse_coast
   coastal_tools.plot_region_box

   depth.add_depth
   depth.add_zmid
   depth.write_time_varying_zmid
   depth.compute_depth
   depth.compute_zmid

   compute_barotropic_streamfunction
   shift_barotropic_streamfunction

.. currentmodule:: mpas_tools.ocean.inject_bathymetry

.. autosummary::
   :toctree: generated/

   inject_bathymetry

.. currentmodule:: mpas_tools.ocean.inject_meshDensity

.. autosummary::
   :toctree: generated/

   inject_meshDensity_from_file
   inject_spherical_meshDensity
   inject_planar_meshDensity

.. currentmodule:: mpas_tools.ocean.inject_preserve_floodplain

.. autosummary::
   :toctree: generated/

   inject_preserve_floodplain


.. currentmodule:: mpas_tools.ocean.viz.transect

.. autosummary::
   :toctree: generated/

   compute_transect
   find_transect_levels_and_weights
   interp_mpas_to_transect_cells
   interp_mpas_to_transect_nodes
   interp_transect_grid_to_transect_nodes
   plot_feature_transects
   plot_transect

.. currentmodule:: mpas_tools.ocean.viz

.. autosummary::
   :toctree: generated/

   add_inset

Sea-ice Tools
=============

.. currentmodule:: mpas_tools.seaice.mask

.. autosummary::
   :toctree: generated/

    extend_seaice_mask

.. currentmodule:: mpas_tools.seaice.mesh

.. autosummary::
   :toctree: generated/

    write_scrip_file
    write_2D_scripfile
    make_mpas_scripfile_on_cells
    make_mpas_scripfile_on_vertices

.. currentmodule:: mpas_tools.seaice.partition

.. autosummary::
   :toctree: generated/

    gen_seaice_mesh_partition
    prepare_partitions
    create_partitions

.. currentmodule:: mpas_tools.seaice.regions

.. autosummary::
   :toctree: generated/

    make_regions_file

.. currentmodule:: mpas_tools.seaice.regrid

.. autosummary::
   :toctree: generated/

    regrid_to_other_mesh

Logging
=======

.. currentmodule:: mpas_tools.logging

.. autosummary::
   :toctree: generated/

   check_call
   LoggingContext

Transects
=========

.. currentmodule:: mpas_tools.transects

.. autosummary::
   :toctree: generated/

   subdivide_great_circle
   cartesian_to_great_circle_distance
   subdivide_planar
   lon_lat_to_cartesian
   cartesian_to_lon_lat

Vector
======

.. currentmodule:: mpas_tools.vector

.. autosummary::
   :toctree: generated/

   Vector
   Vector.angular_distance
   Vector.intersects
   Vector.intersection
   Vector.straddles
   Vector.dot
   Vector.cross
   Vector.det
   Vector.mag

Visualization
=============

.. currentmodule:: mpas_tools.viz.paraview_extractor

.. autosummary::
   :toctree: generated/

   extract_vtk

.. currentmodule:: mpas_tools.viz.mesh_to_triangles

.. autosummary::
   :toctree: generated/

   mesh_to_triangles

.. currentmodule:: mpas_tools.viz.transect

.. autosummary::
   :toctree: generated/

   find_planar_transect_cells_and_weights
   find_spherical_transect_cells_and_weights
   interp_mpas_horiz_to_transect_nodes
   make_triangle_tree
   mesh_to_triangles

.. currentmodule:: mpas_tools.viz.colormaps

.. autosummary::
   :toctree: generated/

   register_sci_viz_colormaps
