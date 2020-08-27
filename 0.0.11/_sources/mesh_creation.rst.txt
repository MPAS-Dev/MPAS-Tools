.. _mesh_creation:

*************
Mesh Creation
*************

Building a Mesh
===============

The :py:func:`mpas_tools.mesh.creation.build_mesh.build_mesh` function is used
create an MPAS mesh using the `JIGSAW <https://github.com/dengwirda/jigsaw>`_
and `JIGSAW-Python <https://github.com/dengwirda/jigsaw-python>`_ packages.

The user must define a local python module ``define_base_mesh`` that provides a
function that returns a 2D array ``cellWidth`` of cell sizes in kilometers.

If the mesh is on a sphere, this function is called ``cellWidthVsLatLon()``
and also returns 1D ``lon`` and ``lat`` arrays.

The mesh is planar, the function is called ``cellWidthVsXY()`` and returns 4
arrays in addition to ``cellWidth``: 1D ``x`` and ``y`` arrays defining planar
coordinates in meters; as well as ``geom_points``, list of point coordinates for
bounding polygon for the planar mesh; and ``geom_edges``, list of edges between
points in ``geom_points`` that define the bounding polygon.

The result is an MPAS mesh file ``base_mesh.nc`` as well as several intermediate
files: ``mesh.log``, ``mesh-HFUN.msh``, ``mesh.jig``, ``mesh-MESH.msh``,
``mesh.msh``, and ``mesh_triangles.nc``.

The :py:func:`mpas_tools.viz.paraview_extractor.extract_vtk` function is used
to produce a VTK file in the ``base_mesh_vtk`` directory that can be viewed in
`ParaVeiw <https://www.paraview.org/>`_.

Optionally, a field, ``cellSeedMask``, can be added to the mesh file that can
later be used preserve a "flood plain" of positive elevations in the MPAS mesh.
See :py:func:`mpas_tools.mesh.creation.inject_preserve_floodplain.inject_preserve_floodplain`.

Optioanlly, a field, ``bottomDepthObserved``, can be added to the mesh file
with bathymetry data from one of two topography files: ``earth_relief_15s.nc``
or ``topo.msh``. If bathymetry should be added to the mesh, a local link with
one of these file names must exist. See
py:func:`mpas_tools.mesh.creation.inject_bathymetry.inject_bathymetry``.

A simple example of ``define_base_mesh.py`` for a spherical mesh with constant,
240-km resolution is:

.. code-block:: python

  import numpy as np


  def cellWidthVsLatLon():
      """
      Create cell width array for this mesh on a regular latitude-longitude grid.

      Returns
      -------
      cellWidth : numpy.ndarray
          m x n array of cell width in km
      lon : numpy.ndarray
          longitude in degrees (length n and between -180 and 180)
      lat : numpy.ndarray
          longitude in degrees (length m and between -90 and 90)
      """

      ddeg = 10
      constantCellWidth = 240

      lat = np.arange(-90, 90.01, ddeg)
      lon = np.arange(-180, 180.01, ddeg)

      cellWidth = constantCellWidth * np.ones((lat.size, lon.size))
      return cellWidth, lon, lat

With this module defined locally, a mesh can be generated either with the
command-line tool ``build_mesh``:

.. code-block::

  $ build_mesh

or by calling py:func:`mpas_tools.mesh.creation.build_mesh.build_mesh`:

.. code-block:: python

  from mpas_tools.mesh.creation.build_mesh import build_mesh


  build_mesh()

The full usage details of the command-line tool are:

.. code-block::

  $ build_mesh --help

  usage: build_mesh [-h] [--preserve_floodplain]
                    [--floodplain_elevation FLOODPLAIN_ELEVATION]
                    [--inject_bathymetry] [--geometry GEOMETRY]
                    [--plot_cellWidth]

  optional arguments:
    -h, --help            show this help message and exit
    --preserve_floodplain
                          Whether a flood plain (bathymetry above z = 0) should
                          be preserved in the mesh
    --floodplain_elevation FLOODPLAIN_ELEVATION
                          The elevation in meters to which the flood plain is
                          preserved, default is 20 m
    --inject_bathymetry   Whether one of the default bathymetry datasets,
                          earth_relief_15s.nc or topo.msh, should be added to
                          the MPAS mesh
    --geometry GEOMETRY   Whether the mesh is on a sphere or a plane, default is
                          a sphere
    --plot_cellWidth      Whether to produce a plot of cellWidth



