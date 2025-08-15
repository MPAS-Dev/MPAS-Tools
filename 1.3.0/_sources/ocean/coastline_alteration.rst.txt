.. _ocean_coastline_alteration:

Coastline Alteration
====================

The :py:mod:`mpas_tools.ocean.coastline_alteration` module contains several
functions for modifying the coastline in MPAS cells.  This is accomplished by
modifying the land mask that indicates which MPAS cells are ocean vs. land
(or grounded ice).

.. _coastline_land_blockages:

Adding Land Blockages
---------------------

The function
:py:func:`mpas_tools.ocean.coastline_alteration.add_critical_land_blockages()`
add the masks associated with one or more transects to the land mask.  This is
useful if there are land features (typically peninsulas or narrow land bridges)
that must block the ocean flow even if they are too narrow to be resolved in
the MPAS mesh.

An example of the typical workflow that uses this function would be:

.. code-block:: python

    import xarray

    from geometric_features import GeometricFeatures
    from mpas_tools.cime.constants import constants
    from mpas_tools.mesh import mask
    from mpas_tools.ocean.coastline_alteration import add_critical_land_blockages


    earthRadius = constants['SHR_CONST_REARTH']

    # an object used to read feature collections from the geometric_features
    # package
    gf = GeometricFeatures()

    # get geometry of the land coverage from Natural Earth
    fcLandCoverage = gf.read(componentName='natural_earth', objectType='region',
                             featureNames=['Land Coverage'])

    # read in the base mesh
    dsBaseMesh = xarray.open_dataset('base_mesh.nc')
    # create a mask on the base mesh that is ones for land or grounded ice and
    # zeros for ocean
    dsLandMask = mask.compute_mpas_region_masks(dsBaseMesh,
                                                fcMask=fcLandCoverage,
                                                maskTypes=('cell',))

    # get a collection of features from geometric_features that are meant for
    # use as critical land blockages
    fcCritBlockages = gf.read(componentName='ocean', objectType='transect',
                              tags=['Critical_Land_Blockage'])

    # make masks from the transects
    dsCritBlockMask = mask.compute_mpas_transect_masks(dsBaseMesh,
                                                       fcMask=fcCritBlockages,
                                                       earthRadius=earthRadius)

    # add the masks to the "land" mask
    dsLandMask = add_critical_land_blockages(dsLandMask, dsCritBlockMask)

.. _coastline_mask_land_locked:

Masking Land-locked Cells
-------------------------

By default, the land mask produced by the MPAS :ref:`mask_creator` can produce
ocean cells with vertices that are all on the land-ocean boundary.  MPAS-Seaice
uses a so-called
`Arakawa B-grid <https://doi.org/10.1016%2FB978-0-12-460817-7.50009-4>`_,
meaning that velocities are located at vertices of the MPAS mesh.  This means
that sea-ice flow into or out of a given cell is not possible if all vertices
of that cell are boundary vertices (where the velocity is, by definition, zero).
This problem is alleviated by calling
:py:func:`mpas_tools.ocean.coastline_alteration.add_land_locked_cells_to_mask()`.
Any "land-locked" cells with only boundary vertices are removed from the mesh.
Land-locked cells are only removed poleward of a threshold latitude (43 degrees
by default).  The user can specify the number of iterations (``nSweeps``) of
land-locked cell removal, since removing land-locked cells can produce new
land-locked cells.

Here is an example workflow that removes land-locked cells:

.. code-block:: python

    import xarray

    from geometric_features import GeometricFeatures
    from mpas_tools.mesh import mask
    from mpas_tools.ocean.coastline_alteration import \
        add_land_locked_cells_to_mask


    # an object used to read feature collections from the geometric_features
    # package
    gf = GeometricFeatures()

    # get geometry of the land coverage from Natural Earth
    fcLandCoverage = gf.read(componentName='natural_earth', objectType='region',
                             featureNames=['Land Coverage'])

    # read in the base mesh
    dsBaseMesh = xarray.open_dataset('base_mesh.nc')
    # create a mask on the base mesh that is ones for land or grounded ice and
    # zeros for ocean
    dsLandMask = mask.compute_mpas_region_masks(dsBaseMesh,
                                                fcMask=fcLandCoverage,
                                                maskTypes=('cell',))

    # Find ocean cells that are land-locked, and alter the cell mask so that
    # they are counted as land cells
    dsLandMask = add_land_locked_cells_to_mask(dsLandMask, dsBaseMesh,
                                               latitude_threshold=43.0,
                                               nSweeps=20)

.. _coastline_widen_transects:

Widening Transects
------------------

Similarly to :ref:`coastline_mask_land_locked`, if critical passages in polar
regions are too narrow, they can become blocked by sea ice that cannot be
advected.  Sea-ice flow is not possible unless channels are at least 2 cells
wide. This widening is accomplished with
:py:func:`mpas_tools.ocean.coastline_alteration.widen_transect_edge_masks()`.
Channels are only widened poleward of a threshold latitude (43 degrees by
default).

An example workflow that includes transect-widening is:

.. code-block:: python

    import xarray

    from geometric_features import GeometricFeatures
    from mpas_tools.cime.constants import constants
    from mpas_tools.mesh import mask
    from mpas_tools.ocean.coastline_alteration import widen_transect_edge_masks

    earthRadius = constants['SHR_CONST_REARTH']

    # an object used to read feature collections from the geometric_features
    # package
    gf = GeometricFeatures()

    # read in the base mesh
    dsBaseMesh = xarray.open_dataset('base_mesh.nc')

    # merge transects for critical passages into critical_passages.geojson
    fcCritPassages = gf.read(componentName='ocean', objectType='transect',
                             tags=['Critical_Passage'])

    # create masks from the transects
    dsCritPassMask = mask.compute_mpas_transect_masks(dsBaseMesh,
                                                      fcMask=fcCritPassages,
                                                      earthRadius=earthRadius)

    # Alter critical passages to be at least two cells wide, to avoid sea ice
    # blockage.
    dsCritPassMask = widen_transect_edge_masks(dsCritPassMask, dsBaseMesh,
                                               latitude_threshold=43.0)

