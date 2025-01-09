.. _ocean_moc:

Meridional Overturning Circulation
==================================

The :py:mod:`mpas_tools.ocean.moc` module contains tools for setting up masks
used by both
`MPAS-Ocean <https://github.com/MPAS-Dev/MPAS-Model/tree/ocean/develop>`_ and
`MPAS-Analysis <https://mpas-dev.github.io/MPAS-Analysis/stable/>`_ for
computing the global and basin-wide Meridional Overturning Circulation (MOC).

.. _moc_basins:

Building MOC Basins
-------------------

The :py:func:`mpas_tools.ocean.moc.build_moc_basins()` function can be used to
merge a large number of individual seas from the
`geometric_features <http://mpas-dev.github.io/geometric_features/stable/>`
package into 5 larger features defining the MOC basins:

* Global
* Atlantic
* Indian
* Pacific
* Indo-pacific

A call to this function returns a
`FeatureCollection <http://mpas-dev.github.io/geometric_features/stable/feature_collection.html>`_
with these 5 basins.

A typical workflow calling this function would look like:

.. code-block:: python

    from geometric_features import GeometricFeatures
    from mpas_tools.ocean.moc import build_moc_basins


    gf = GeometricFeatures()
    fcMOC = build_moc_basins(gf)

.. _moc_southern_transects:

Adding Southern Transects
-------------------------

Typically, the basins on their own are not particularly useful.  For all but
the global MOC, a southern transect is needed to determine the flow into or
out of the basin from the south.  The function
:py:func:`mpas_tools.ocean.mocadd_moc_southern_boundary_transects()` can be
used to add masks for all the edges associated with the southern boundary of
each MOC basin.

A typical workflow calling this function would look like:

.. code-block:: python

    import xarray
    from geometric_features import GeometricFeatures
    from mpas_tools.ocean.moc import build_moc_basins


    mesh_filename = 'mesh.nc'

    gf = GeometricFeatures()
    fcMOC = build_moc_basins(gf)

    dsMesh = xarray.open_dataset(mesh_filename)
    dsMasks = mpas_tools.mesh.conversion.mask(dsMesh=dsMesh, fcMask=fcMOC)

    dsMasksAndTransects = add_moc_southern_boundary_transects(dsMasks, dsMesh)

In this example, only the ``mesh.nc`` file is required as an input.  The
resulting ``xarray.Dataset`` contains both the basin and southern-transect
masks.

A command-line tool ``moc_southern_boundary_extractor`` is also available
for this purpose:

.. code-block:: none

    $ moc_southern_boundary_extractor --help

    usage: moc_southern_boundary_extractor [-h] -f IN_FILE -m MESH_FILE -o
                                              OUT_FILE

    This script takes a mesh file (-m flag) and a file with MOC regions masks
    (-f flag) produce by the MPAS mask creator.  The script produces a copy of
    the contents of the MOC mask file, adding transects that mark the southern
    boundary of each region in a file indicated with the -o flag.  The transect
    is applied only to vertices and edges, not cells, because the need for southern
    boundary transect data on cells is not foreseen.

    optional arguments:
      -h, --help            show this help message and exit
      -f IN_FILE, --in_file IN_FILE
                            Input file with MOC masks
      -m MESH_FILE, --mesh_file MESH_FILE
                            Input mesh file
      -o OUT_FILE, --out_file OUT_FILE
                            Output file for MOC masks and southern-boundary transects

The command-line tool is largely intended for backwards compatibility and the
python function is the preferred way of building a workflow with this
functionality.

.. _moc_basins_and_transects:

Building MOC Basins and Transects Together
------------------------------------------

Typically, a workflow can be made more efficient by using the function
:py:func:`mpas_tools.ocean.moc.make_moc_basins_and_transects()` takes care of
both :ref:`moc_basins` and :ref:`moc_southern_transects` as well as writing out
the results to files.

A typical workflow calling this function would look like:

.. code-block:: python

    from geometric_features import GeometricFeatures
    from mpas_tools.ocean.moc import make_moc_basins_and_transects


    mesh_filename = 'mesh.nc'
    mesh_name = 'EC30to60kmL60E3SMv2r03'

    mask_filename = '{}_moc_masks.nc'.format(mesh_name)
    mask_and_transect_filename = '{}_moc_masks_and_transects.nc'.format(
        mesh_name)

    geojson_filename = 'moc_basins.geojson'

    gf = GeometricFeatures()

    make_moc_basins_and_transects(gf, mesh_filename, mask_and_transect_filename,
                                  geojson_filename=geojson_filename,
                                  mask_filename=mask_filename)

In this example, only the ``mesh.nc`` file is required as an input.  The basin
and transect masks are written to ``mask_and_transect_filename``, and we also
request that the intermediate data sets get written out for, perhaps for
purposes of debugging or provenance.  The MOC feature collection from
:ref:`moc_basins` will be written to ``geojson_filename``, while the basin
masks (without the associate transect masks) will be written to
``mask_filename``.
