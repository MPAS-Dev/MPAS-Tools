.. _mesh_creation:

*************
Mesh Creation
*************

Building a Mesh
===============

The :py:mod:`mpas_tools.mesh.creation.build_mesh` module is used
create an MPAS mesh using the `JIGSAW <https://github.com/dengwirda/jigsaw>`_
and `JIGSAW-Python (jigsawpy) <https://github.com/dengwirda/jigsaw-python>`_
packages.


Spherical Meshes
----------------

Spherical meshes are constructed with the function
:py:func:`mpas_tools.mesh.creation.build_mesh.build_spherical_mesh()`.
The user provides a 2D array ``cellWidth`` of cell sizes in kilometers along
1D arrays for the longitude and latitude (the cell widths must be on a lon/lat
tensor grid) and the radius of the earth in meters.

The result is an MPAS mesh file, called ``base_mesh.nc`` by default, as well as
several intermediate files: ``mesh.log``, ``mesh-HFUN.msh``, ``mesh.jig``,
``mesh-MESH.msh``, ``mesh.msh``, and ``mesh_triangles.nc``.

Here is a simple example script for creating a uniform MPAS mesh with 240-km
resolution:

.. code-block:: python

    #!/usr/bin/env python
    import numpy as np
    from mpas_tools.ocean import build_spherical_mesh


    def cellWidthVsLatLon():
        """
        Create cell width array for this mesh on a regular latitude-longitude grid.
        Returns
        -------
        cellWidth : ndarray
            m x n array of cell width in km
        lon : ndarray
            longitude in degrees (length n and between -180 and 180)
        lat : ndarray
            longitude in degrees (length m and between -90 and 90)
        """
        dlat = 10
        dlon = 10
        constantCellWidth = 240

        nlat = int(180/dlat) + 1
        nlon = int(360/dlon) + 1

        lat = np.linspace(-90., 90., nlat)
        lon = np.linspace(-180., 180., nlon)

        cellWidth = constantCellWidth * np.ones((lat.size, lon.size))
        return cellWidth, lon, lat


    def main():
        cellWidth, lon, lat = cellWidthVsLatLon()
        build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc')


    if __name__ == '__main__':
        main()

We define the resolution on a coarse (10 degree by 10 degree) grid because it
is uniform.  Meshes with more complex variation may require higher resolution
grids to cell widths.

Planar Meshes
-------------

Planar meshes can be constructed with the function
:py:func:`mpas_tools.mesh.creation.build_mesh.build_planar_mesh()`.  Provide
this function with a 2D array ``cellWidth`` of cell sizes in kilometers and
1D arrays for x and y (the cell widths must be on a 2D tensor grid).  Planar
meshes also require ``geom_points``, a list of point coordinates for bounding
polygon for the planar mesh, and ``geom_edges``, a list of edges between points
in ``geom_points`` that define the bounding polygon.

As for spehrical meshes, the result is an MPAS mesh file, called
``base_mesh.nc`` by default, as well as several intermediate files:
``mesh.log``, ``mesh-HFUN.msh``, ``mesh.jig``, ``mesh-MESH.msh``, ``mesh.msh``,
and ``mesh_triangles.nc``.


JIGSAW Driver
-------------

Underlying both spherical and planar mesh creation is the JIGSAW driver
function :py:func:`mpas_tools.mesh.creation.jigsaw_driver.jigsaw_driver()`.  This
function is used to setup data structures and then build a JIGSAW mesh using
``jigsawpy``.


Converting Between Mesh Formats
===============================

MSH to MPAS NetCDF
------------------

``jigsawpy`` produces meshes in ``.msh`` format that need to be converted to
`NetCDF <https://www.unidata.ucar.edu/software/netcdf/>`_ files for use by MPAS
components.  A utility function
:py:func:`mpas_tools.mesh.creation.jigsaw_to_netcdf.jigsaw_to_netcdf()` or the
command-line utility ``jigsaw_to_netcdf`` are used for this purpose.

In addition to the input ``.msh`` and output ``.nc`` files, the user must
specify whether this is a spherical or planar mesh and, if it is spherical,
provide the radius of the Earth in meters.

Triangle to MPAS NetCDF
-----------------------

Meshes in `Triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_ format
can be converted to MPAS NetCDF format using
:py:func:`mpas_tools.mesh.creation.triangle_to_netcdf.triangle_to_netcdf()` or
the ``triangle_to_netcdf`` command-line tool.

The user supplies the names of input ``.node`` and ``.ele`` files and the
name of an output MPAS mesh file.

MPAS NetCDF to Triangle
-----------------------

MPAS meshes in NetCDF format can be converted to ``Triangle`` format using
:py:func:`mpas_tools.mesh.creation.mpas_to_triangle.mpas_to_triangle()` or
the ``mpas_to_triangle`` command-line tool.

The user supplies the name of an input MPAS mesh file and the output prefix
for the resulting Triangle ``.node`` and ``.ele`` files.

MPAS NetCDF to SCRIP
--------------------

The function :py:func:`mpas_tools.scrip.from_mpas.scrip_from_mpas()` can be
used to convert an MPAS mesh file in NetCDF format to
`SCRIP <http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_5_2_0rp1/ESMF_refdoc/node3.html#SECTION03024000000000000000>`_
format.  SCRIP files are typically used to create mapping files used to
interpolate between meshes.
