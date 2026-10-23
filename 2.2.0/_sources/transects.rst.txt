.. _transects:

*********
Transects
*********

The :py:mod:`mpas_tools.transects` module contains functions used to define
transects through MPAS meshes.  These transects can be used to create masks
of the cells, edges or dual-mesh cells (in some sense "vertices") that
intersect the transect.  They can also be used for visualization such as
plotting vertical cross-sections of MPAS data along the transect.

.. _subdividing_transects:

Subdividing transects
=====================

For both visualization and intersection detection, it is often useful to
subdivide a transect into smaller segments.  This is performed with the
function :py:func:`mpas_tools.transects.subdivide_great_circle()` for spherical
meshes and with :py:func:`mpas_tools.transects.subdivide_planar()` for planar
meshes.

For spherical meshes, subdivision is performed in Cartesian coordinates.  Since
transects are typically provided as a sequence of longitude/latitude points,
it is typically necessary to convert to Cartesian coordinates using
:py:func:`mpas_tools.transects.lon_lat_to_cartesian()` and then back to
longitude/latitude coordinates using
:py:func:`mpas_tools.transects.cartesian_to_lon_lat()`.


Low-level functions
===================

The module also shares some lower-level functions used elsewhere in the
package.

The arc length (in radians) along a transect can be found with
:py:func:`mpas_tools.transects.angular_distance()`.

The function :py:func:`mpas_tools.transects.intersects()` can be used to
determine if 2 arcs on the sphere intersect one another and
:py:func:`mpas_tools.transects.intersection()` can be used to find the
intersection point of 2 intersecting arcs.
