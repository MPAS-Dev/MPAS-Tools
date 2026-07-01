.. _ocean_mesh_creation:

Ocean Mesh Creation
===================

The :py:mod:`mpas_tools.ocean.build_mesh` module is used create
ocean-specific MPAS meshes using the
:py:mod:`mpas_tools.mesh.creation.build_mesh` module.

Spherical meshes are constructed with the function
:py:func:`mpas_tools.ocean.build_mesh.build_spherical_mesh()`.  The basic
arguments are the same as those to
:py:func:`mpas_tools.mesh.creation.build_mesh.build_spherical_mesh()`.

Similarly, planar meshes can be constructed with the function
:py:func:`mpas_tools.ocean.build_mesh.build_planar_mesh()`, which has the
same basic arguments as
:py:func:`mpas_tools.mesh.creation.build_mesh.build_planar_mesh()`.

Each of these functions has additional, optional arguments that allow users to:

  * specify a directory for extracting VTK geometry for viewing in
    `ParaVeiw <https://www.paraview.org/>`_.
    The :py:func:`mpas_tools.viz.paraview_extractor.extract_vtk` function is
    used to produce a VTK file in this directory, named ``base_mesh_vtk``
    by default.

  * Specify whether to preserve a region of the mesh above sea level as a
    floodplain, and the elevation up to which this regions should remain part
    of the mesh.  This feature is used in coastal simulations that support
    wetting and drying.  A field, ``cellSeedMask``, is added to the mesh file
    that can later be used preserve the floodplain.
    See :py:func:`mpas_tools.ocean.inject_preserve_floodplain.inject_preserve_floodplain`.


  * Whether to add a default bathymetry data set to the mesh.  A field,
    ``bottomDepthObserved``, is added to the mesh file with bathymetry data
    from one of two topography files: ``earth_relief_15s.nc`` or ``topo.msh``.
    If bathymetry should be added to the mesh, a local link with one of these
    file names must exist.
    See :py:func:`mpas_tools.ocean.inject_bathymetry.inject_bathymetry`.
