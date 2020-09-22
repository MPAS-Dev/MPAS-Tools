.. _mesh_creation:

.. |---| unicode:: U+2014  .. em dash, trimming surrounding whitespace
   :trim:

*************
Mesh Creation
*************

Uniform, Planar Meshes
======================

The most basic tool for creating an MPAS mesh is the function
:py:func:`mpas_tools.planar_hex.make_planar_hex_mesh()` or the associated
command-line tool ``planar_hex``.  These tools create a uniform, planar mesh
of hexagons that may optionally be periodic in either or both of the x and
y directions.  A typical call to ``make_planar_hex_mesh()`` might look like
the following:

.. code-block:: python

    from mpas_tools.planar_hex import make_planar_hex_mesh

    dsMesh = make_planar_hex_mesh(nx=10, ny=20, dc=1000., nonperiodic_x=False,
                                  nonperiodic_y=False)

This creates a periodic mesh in both x and y that is 10 grid cells by 20 grid
cells in size with 1-km resolution.  The mesh is not save to a file in this
example but is instead returned as an ``xarray.Dataset`` object.

The command-line approach for generating the same mesh would be:

.. code-block:: none

    $ planar_hex --nx=10 --ny=20 --dc=1000. \
         --outFileName='periodic_mesh_10x20_1km.nc

In this case, the resulting mesh is written to a file.

Building a JIGSAW Mesh
======================

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

Mesh Definition Tools
=====================

The :py:mod:`mpas_tools.mesh.creation.mesh_definition_tools` module includes
several tools for defining the ``cellWidth`` variable.

Merging Cell Widths
-------------------
The function
:py:func:`mpas_tools.mesh.creation.mesh_definition_tools.mergeCellWidthVsLat()`
is used to combine two cell-width distributions that are functions of latitude
only and which asymptote to different constant values north and south of a given
transition latitude with a ``tanh`` function of a given characteristic width.

For example, the following code snippet will produce cell widths as a function
of latitude of about 30 km south of the Arctic Circle and 10 km north of that
latitude, transitioning over a characteristic "distance" of about 5 degrees.

.. code-block:: python

    import numpy
    from mpas_tools.mesh.creation.mesh_definition_tools import \
        mergeCellWidthVsLat


    lat = numpy.linspace(-90., 90., 181)
    cellWidthInSouth = 30.
    cellWidthInNorth = 10.
    latTransition = 66.5
    latWidthTransition = 5.

    cellWidths = mergeCellWidthVsLat(lat, cellWidthInSouth, cellWidthInNorth,
        latTransition, latWidthTransition)

.. _ec_mesh:

Defining an Eddy-closure Mesh
-----------------------------

One of the commonly used flavor of MPAS-Ocean and MPAS-Seaice meshes is designed
with relatively coarse resolution in mind (requiring parameterization of ocean
eddies with an "eddy closure").  This flavor of mesh has resolution that is
purely a function of latitude, with 5 regions of relatively uniform resolution
(north polar, northern mid-latitudes, equatorial, southern mid-latitudes and
south polar) with smooth (``tanh``) transitions between these resolutions.

The default EC mesh has resolutions of 35 km at the poles, 60 km at
mid-latitudes and 30 km at the equator.  Transitions between equatorial and
mid-latitude regions are at 15 degrees N/S latitude and transitions between
mid-latitude and polar regions are at 73 degrees N/S latitude.  The
transition near the equator is somewhat more abrupt (~6 degrees) than near the
poles (~9 degrees).  The switch between the transitional ``tanh`` functions is
made at 40 degrees N/S latitude, where the resolution is nearly constant and no
appreciable discontinuity arises.  The default EC mesh can be obtained with the
function
:py:func:`mpas_tools.mesh.creation.mesh_definition_tools.EC_CellWidthVsLat()`:

.. code-block:: python

    import numpy
    from mpas_tools.mesh.creation.mesh_definition_tools import \
        EC_CellWidthVsLat

    lat = numpy.linspace(-90., 90., 181)
    cellWidths = EC_CellWidthVsLat(lat)

.. _rrs_mesh:

Defining a Rossby-radius Mesh
-----------------------------

Another common flavor of MPAS-Ocean and MPAS-Seaice meshes is designed for
higher resolutions, where the Rossby radius of deformation can be (at least
partially) resolved.  These meshes approximately scale their resolution in
proportion to the Rossby radius.

A typical Rossby Radius Scaling (RRS) mesh has a resolution at the poles that is
three times finer than the resolution at the equator.  For example, the RRS mesh
used in E3SMv1 high resolution simulations would be defined, using the function
:py:func:`mpas_tools.mesh.creation.mesh_definition_tools.RRS_CellWidthVsLat()`
by:

.. code-block:: python

    import numpy
    from mpas_tools.mesh.creation.mesh_definition_tools import \
        RRS_CellWidthVsLat

    lat = numpy.linspace(-90., 90., 181)
    cellWidths = RRS_CellWidthVsLat(lat, cellWidthEq=18., cellWidthPole=6.)

Defining an Atlantic/Pacific Mesh
---------------------------------

The function
:py:func:`mpas_tools.mesh.creation.mesh_definition_tools.AtlanticPacificGrid()`
can be used to define a mesh that has two different, constant resolutions in the
Atlantic and Pacific Oceans.


Signed Distance Functions
=========================

The :py:mod:`mpas_tools.mesh.creation.signed_distance` module includes several
functions for creating ``cellWidth`` variables based on the signed distance from
a boundary curve on the sphere.  A signed distance function is positive outside
the bounding shape and negative inside, with a value proportional to the
distance to the nearest point on the curve (so the function is equal to zero on
the curve).  Signed distance functions provide a useful way ot define
transitions in resolution based on complex shapes that can be defined using
`geojson <https://geojson.org/>`_ files.  These files can be created by hand,
e.g. at `geojson.io <http://geojson.io/>`_ or in python using libraries like
`shapely <https://shapely.readthedocs.io/en/stable/index.html>`_.

Calls to the functions in this module require a
`FeatureCollection <http://mpas-dev.github.io/geometric_features/stable/feature_collection.html>`_
object from the
`geometric_features <http://mpas-dev.github.io/geometric_features/stable/index.html>`_
package.  The ``FeatureColleciton`` must define one or more regions on the
sphere from which the distance, mask, or signed distance will be computed.
The ``FeatureColleciton`` could come from the predefined features included in
``geometric_features``, could be read in from a ``geojson`` file (see
`Reading in Features <http://mpas-dev.github.io/geometric_features/stable/feature_collection.html#reading-in-features>`_),
or could be created as part of a python script with ``shapely`` or other tools.

In this example, we first define a base resolution using the default EC mesh
(see :ref:`ec_mesh`) and then use
:py:func:`mpas_tools.mesh.creation.signed_distance.signed_distance_from_geojson()`
to create a signed distance function from a ``FeatureCollection`` read in from
`this geojson file <https://github.com/MPAS-Dev/MPAS-Model/blob/ocean/develop/testing_and_setup/compass/ocean/global_ocean/SO60to10wISC/init/high_res_region.geojson>`_.
The signed distance function is used to define a region of high resolution (12
km) around Antarctica.

.. code-block:: python

    import numpy as np
    import mpas_tools.mesh.creation.mesh_definition_tools as mdt
    from mpas_tools.mesh.creation.signed_distance import \
        signed_distance_from_geojson
    from geometric_features import read_feature_collection
    from mpas_tools.cime.constants import constants


    dlon = 0.1
    dlat = dlon
    earth_radius = constants['SHR_CONST_REARTH']
    nlon = int(360./dlon) + 1
    nlat = int(180./dlat) + 1
    lon = np.linspace(-180., 180., nlon)
    lat = np.linspace(-90., 90., nlat)

    cellWidth = mdt.EC_CellWidthVsLat(lat)

    # broadcast cellWidth to 2D
    _, cellWidth = np.meshgrid(lon, cellWidthVsLat)

    # now, add the high-res region
    fc = read_feature_collection('high_res_region.geojson')

    so_signed_distance = signed_distance_from_geojson(fc, lon, lat,
                                                      earth_radius,
                                                      max_length=0.25)

    # Equivalent to 20 degrees latitude
    trans_width = 1600e3
    trans_start = -500e3
    dx_min = 12.

    weights = 0.5 * (1 + np.tanh((so_signed_distance - trans_start) /
                                 trans_width))

    cellWidth = dx_min * (1 - weights) + cellWidth * weights

Sometimes it can be useful to extract just the mask of the region of interest
(defined as ``0`` outside the the region and ``1`` inside it) or the unsigned
distance.  For these purposes, use the functions
:py:func:`mpas_tools.mesh.creation.signed_distance.mask_from_geojson()`
and
:py:func:`mpas_tools.mesh.creation.signed_distance.distance_from_geojson()`,
respectively.

