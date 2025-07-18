.. _mesh_conversion:

***************
Mesh Conversion
***************

.. _mesh_converter:

Mesh Converter
==============

The ``MpasMeshConverter.x`` command-line tool and the Python wrapper
:py:func:`mpas_tools.mesh.conversion.convert` convert a dataset describing
cell and vertex locations and connectivity into a valid MPAS mesh that
follows the `MPAS mesh specification
<https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf>`_.

Example usage (command line):

.. code-block::

  $ planar_hex --nx 4 --ny 4 --dc 10e3 -o base_mesh.nc
  $ MpasMeshConverter.x base_mesh.nc mesh.nc

This generates a small, doubly periodic MPAS mesh and converts it to a
spec-compliant format. ``MpasMeshConverter.x`` takes input and output mesh
filenames as arguments; if omitted, it prompts for them.

Equivalent Python usage:

.. code-block:: python

  from mpas_tools.planar_hex import make_planar_hex_mesh
  from mpas_tools.mesh.conversion import convert
  from mpas_tools.io import write_netcdf

  ds = make_planar_hex_mesh(nx=4, ny=4, dc=10e3, nonperiodic_x=False,
                            nonperiodic_y=False)
  ds = convert(ds)
  write_netcdf(ds, 'mesh.nc')

**Input requirements:** The mesh must define the following dimensions,
variables, and global attributes (example sizes shown):

.. code-block::

  netcdf mesh {
  dimensions:
    nCells = 16 ;
    nVertices = 32 ;
    vertexDegree = 3 ;
  variables:
    double xCell(nCells) ;
    double yCell(nCells) ;
    double zCell(nCells) ;
    double xVertex(nVertices) ;
    double yVertex(nVertices) ;
    double zVertex(nVertices) ;
    int cellsOnVertex(nVertices, vertexDegree) ;
    double meshDensity(nCells) ;

  // global attributes:
        :on_a_sphere = "NO" ;
        :sphere_radius = 0. ;
        :is_periodic = "YES" ;

The ``meshDensity`` variable is required for historical reasons and is passed
unchanged to the output mesh.

Optional global attributes (passed through):

.. code-block::

  // global attributes:
        :x_period = 40000. ;
        :y_period = 34641.0161513775 ;
        :history = "Tue May 26 20:58:10 2020: /home/xylar/miniconda3/envs/mpas/bin/planar_hex --nx 4 --ny 4 --dc 10e3 -o base_mesh.nc" ;

If present, the ``file_id`` attribute is preserved as ``parent_id`` in the
output mesh, and a new ``file_id`` is generated.

The converter also generates a ``graph.info`` file for graph partitioning
tools (e.g., Metis). In Python, this file is only written if the
``graphInfoFileName`` argument is provided.

.. _cell_culler:

Cell Culler
===========

The ``MpasCellCuller.x`` command-line tool and the Python wrapper
:py:func:`mpas_tools.mesh.conversion.cull` remove cells from a mesh based on
the ``cullCell`` field and/or provided mask datasets. The culling logic is:

- The ``cullCell`` field, mask(s) from a masking dataset, and the inverse of
  mask(s) from an inverse-masking dataset are merged (union).
- A preserve-masking dataset indicates cells that must *not* be culled.

Example workflow (command line):

.. code-block::

  $ merge_features -c natural_earth -b region -n "Land Coverage" -o land.geojson
  $ MpasMaskCreator.x base_mesh.nc land.nc -f land.geojson
  $ MpasCellCuller.x base_mesh.nc culled_mesh.nc -m land.nc

This merges features to create a land mask, generates a mask on the mesh,
and culls cells where the mask is 1.

Equivalent Python workflow:

.. code-block:: python

  import xarray
  from geometric_features import GeometricFeatures
  from mpas_tools.mesh.conversion import mask, cull

  gf = GeometricFeatures()
  fcLandCoverage = gf.read(
      componentName='natural_earth',
      objectType='region',
      featureNames=['Land Coverage']
  )
  dsBaseMesh = xarray.open_dataset('base_mesh.nc')
  dsLandMask = mask(dsBaseMesh, fcMask=fcLandCoverage)
  dsCulledMesh = cull(dsBaseMesh, dsMask=dsLandMask)
  write_netcdf(dsCulledMesh, 'culled_mesh.nc')

Full usage of ``MpasCellCuller.x``:

.. code-block::

    MpasCellCuller.x [input_name] [output_name] [[-m/-i/-p] masks_name] [-c]

        input_name:         Input MPAS mesh.
        output_name:        Output culled MPAS mesh (default: culled_mesh.nc).
        -m/-i/-p:           Masking options:
            -m: Mask file(s) (1 = cull cell).
            -i: Inverse mask file(s) (0 = cull cell).
            -p: Preserve mask file(s) (1 = do not cull cell).
        -c:                 Output cell mapping files.

.. _mask_creator:

Mask Creator
============

The ``MpasMaskCreator.x`` command-line tool and the Python wrapper
:py:func:`mpas_tools.mesh.conversion.mask` create region masks from features
or seed points.

Example usage is shown above under Cell Culler.

Full usage of ``MpasMaskCreator.x``:

.. code-block::

    MpasMaskCreator.x in_file out_file [ [-f/-s] file.geojson ] [--positive_lon]
        in_file: Input mesh file.
        out_file: Output mask file.
        -s file.geojson: Use points as seed locations for flood fill.
        -f file.geojson: Use features (regions, transects, or points) for masks.
        --positive_lon: Use 0-360 longitude range for non-standard geojson files.

.. note::
    Temporary files are created and deleted automatically by the Python wrappers.
    Command-line tools require the relevant executables to be available in the path.

.. _py_mask_creation:

Mask Creation with Python Multiprocessing
=========================================

The ``mpas_tools.mesh.mask`` module provides a set of Python functions for
creating region and transect masks on MPAS meshes and longitude/latitude grids.
These functions are designed to be more efficient and flexible than the legacy
serial C++ Mask Creator, especially when used with Python's multiprocessing.

Key Functions
-------------

+-----------------------------------------------+-------------------------------------------------------------+
| Function                                      | Purpose                                                     |
+===============================================+=============================================================+
| compute_mpas_region_masks                     | Create region masks (polygons) on MPAS meshes               |
+-----------------------------------------------+-------------------------------------------------------------+
| compute_mpas_transect_masks                   | Create transect masks (lines) on MPAS meshes                |
+-----------------------------------------------+-------------------------------------------------------------+
| compute_mpas_flood_fill_mask                  | Create a mask by flood-filling from seed points             |
+-----------------------------------------------+-------------------------------------------------------------+
| compute_lon_lat_region_masks                  | Create region masks on a 2D lon/lat grid                    |
+-----------------------------------------------+-------------------------------------------------------------+
| compute_projection_grid_region_masks          | Create region masks on a projected (e.g., polar) grid       |
+-----------------------------------------------+-------------------------------------------------------------+

All of these functions accept a ``pool`` argument (a ``multiprocessing.Pool``)
to parallelize the computation, which is highly recommended for large meshes or
grids. If ``pool=None``, the computation will be performed serially, which may
be slow for large datasets.

General Usage
-------------

The typical workflow is:

1. Open your MPAS mesh or grid as an ``xarray.Dataset``.
2. Read a ``geometric_features.FeatureCollection`` (e.g., from a GeoJSON file).
3. Optionally, create a multiprocessing pool using
   :py:func:`mpas_tools.parallel.create_pool`.
4. Call the appropriate mask creation function, passing the mesh/grid, feature
   collection, and pool.
5. Write the resulting masks to a NetCDF file using
   :py:func:`mpas_tools.io.write_netcdf`.

Example: Creating Region Masks on an MPAS Mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import xarray as xr
    from geometric_features import read_feature_collection
    from mpas_tools.mesh.mask import compute_mpas_region_masks
    from mpas_tools.parallel import create_pool
    from mpas_tools.io import write_netcdf

    dsMesh = xr.open_dataset('mesh.nc', decode_cf=False, decode_times=False)
    fcMask = read_feature_collection('regions.geojson')
    pool = create_pool(process_count=8)
    dsMasks = compute_mpas_region_masks(
        dsMesh, fcMask, maskTypes=('cell', 'vertex'), pool=pool
    )
    write_netcdf(dsMasks, 'region_masks.nc')

Arguments and Options
---------------------

All mask creation functions share several common arguments:

- ``logger``: Optional logger for progress output.
- ``pool``: Optional multiprocessing pool for parallel computation.
- ``chunkSize``: Number of points to process per chunk (default: 1000).
- ``showProgress``: Whether to display a progress bar.
- ``subdivisionThreshold`` or ``subdivisionResolution``: Controls subdivision
  of large polygons or transects for efficiency.

Refer to the Python docstrings or the command-line ``--help`` output for
details on each function's arguments.

Performance Note
----------------

For large meshes or grids, using a multiprocessing pool (via the ``pool``
argument) is strongly recommended for reasonable performance. The pool should
be created early in your script, before large objects are loaded into memory,
and terminated when no longer needed.

Extensibility and Limitations
-----------------------------

- The masking functions are extensible and can be adapted for new types of
  features or grids.
- The algorithms use the ``shapely`` library for geometric operations, which
  is designed for 2D Cartesian geometry. Care is taken to handle longitude
  periodicity, but there may be limitations near the poles or for very large
  polygons.
- For advanced use cases (e.g., custom mask types or additional properties),
  see the source code and docstrings for guidance.

See also the API documentation for :py:mod:`mpas_tools.mesh.mask` for further details.

See also the API documentation for :py:mod:`mpas_tools.mesh.mask` for further details.
                            An MPAS mesh file
      -g GEOJSON_FILE_NAME, --geojson_file_name GEOJSON_FILE_NAME
                            An Geojson file containing mask regions
      -o MASK_FILE_NAME, --mask_file_name MASK_FILE_NAME
                            An output MPAS region masks file
      -t MASK_TYPES [MASK_TYPES ...], --mask_types MASK_TYPES [MASK_TYPES ...]
                            Which type(s) of masks to make: cell, edge or vertex.
                            Default is cell and vertex.
      -c CHUNK_SIZE, --chunk_size CHUNK_SIZE
                            The number of cells, vertices or edges that are
                            processed in one operation
      --show_progress       Whether to show a progress bar
      -s SUBDIVISION, --subdivision SUBDIVISION
                            A threshold in degrees (lon or lat) above which the
                            mask region will be subdivided into smaller polygons
                            for faster intersection checking
      --process_count PROCESS_COUNT
                            The number of processes to use to compute masks. The
                            default is to use all available cores
      --multiprocessing_method MULTIPROCESSING_METHOD
                            The multiprocessing method use for python mask
                            creation ('fork', 'spawn' or 'forkserver')


Computing Transect Masks
------------------------

The function :py:func:`mpas_tools.mesh.mask.compute_mpas_transect_masks()`
and the ``compute_mpas_transect_masks`` command-line tool
are similar to the function for computing region masks.  The function takes a
:py:class:`geometric_features.FeatureCollection` ``fcMask`` that is made up of
transects, rather than regions.  One mask is produced for each feature in the
collection, indicating where the transect
intersects the cell, edge or vertex polygons (see the
`MPAS Mesh Specification <https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf>`_).

The arguments ``logger``, ``pool``, ``chunkSize`` and ``showProgress`` are the
same as for region-mask creation above.

The argument ``subdivisionResolution`` is a length in meters, above which
segments of the transect are subdivided to provide a better representation of
the spherical path in longitude/latitude space.  The default value of 10 km is
typically good enough to capture distortion at typical MPAS mesh resolutions.

The algorithm perform intersections in longitude/latitude space using the
``shapely`` library.  Because ``shapely`` is designed for 2D shapes in a
Cartesian plane, it is not designed for spherical coordinates.  Care has been
taken to handle periodicity at the dateline (antimeridian) but there may be
issues with MPAS mesh polygons containing the north or south pole.  If a user
needs to handle a transect that is very close to the pole, it is likely worth
contacting the developers to request modifications to the code to support this
case.

The resulting variables are:

  - ``transectCellMasks(nCells, nTransects)`` - a cell mask (1 if the transect
    intersects the cell and 0 if not) for each transect
  - ``transectEdgeMasks(nEdges, nTransects)`` - an edge mask for each transect
  - ``transectVertexMasks(nVertices, nTransects)`` - a vertex mask for each
    transect
  - ``transectNames(nTransects, string64)`` - the names of the transects

We don't currently provide cell, edge or vertex indices (e.g.
``transectCellGlobalIDs``) for path along a transect.  This is, in part,
because the algorithm doesn't keep track of the relative order of points along
a transect. This could be updated in the future if there is sufficient demand.

The edge sign (``transectEdgeMaskSigns``) is computed only if
``addEdgeSign=True``, since this takes extra time to compute and isn't always
needed.

.. note::

    While the default ``subdivisionResolution`` is 10 km for
    :py:func:`mpas_tools.mesh.mask.compute_mpas_transect_masks()`, the default
    behavior in the command-line tool ``compute_mpas_transect_masks`` is no
    subdivision because there is otherwise not a good way to specify at the
    command line that no subdivision is desired.  Typically, users will want
    to request subdivision with something like ``-s 10e3``

The command-line tool takes the following arguments:

.. code-block::

    $ compute_mpas_transect_masks --help
    usage: compute_mpas_transect_masks [-h] -m MESH_FILE_NAME -g GEOJSON_FILE_NAME
                                       -o MASK_FILE_NAME
                                       [-t MASK_TYPES [MASK_TYPES ...]]
                                       [-c CHUNK_SIZE] [--show_progress]
                                       [-s SUBDIVISION]
                                       [--process_count PROCESS_COUNT]
                                       [--multiprocessing_method MULTIPROCESSING_METHOD]

    optional arguments:
      -h, --help            show this help message and exit
      -m MESH_FILE_NAME, --mesh_file_name MESH_FILE_NAME
                            An MPAS mesh file
      -g GEOJSON_FILE_NAME, --geojson_file_name GEOJSON_FILE_NAME
                            An Geojson file containing transects
      -o MASK_FILE_NAME, --mask_file_name MASK_FILE_NAME
                            An output MPAS transect masks file
      -t MASK_TYPES [MASK_TYPES ...], --mask_types MASK_TYPES [MASK_TYPES ...]
                            Which type(s) of masks to make: cell, edge or vertex.
                            Default is cell, edge and vertex.
      -c CHUNK_SIZE, --chunk_size CHUNK_SIZE
                            The number of cells, vertices or edges that are
                            processed in one operation
      --show_progress       Whether to show a progress bar
      -s SUBDIVISION, --subdivision SUBDIVISION
                            The maximum resolution (in meters) of segments in a
                            transect. If a transect is too coarse, it will be
                            subdivided. Default is no subdivision.
      --process_count PROCESS_COUNT
                            The number of processes to use to compute masks. The
                            default is to use all available cores
      --multiprocessing_method MULTIPROCESSING_METHOD
                            The multiprocessing method use for python mask
                            creation ('fork', 'spawn' or 'forkserver')
      --add_edge_sign       Whether to add the transectEdgeMaskSigns variable


Computing a Flood-fill Mask
---------------------------

The function :py:func:`mpas_tools.mesh.mask.compute_mpas_flood_fill_mask()`
and the command-line tool ``compute_mpas_flood_fill_mask``
fill in a mask, starting with the cell centers closest to the seed points
given in :py:class:`geometric_features.FeatureCollection` ``fcSeed``.  This
algorithm runs in serial, and will be more efficient the more seed points
are provided and the more widely scattered over the mesh they are.

An optional ``daGrow`` argument to the function (not currently available from
the command-line tool) provides a mask into which the flood fill is allowed to
grow.  The default is all ones.

The resulting dataset contains a single variable:

  - ``cellSeedMask(nCells)`` - a cell mask that is 1 where the flood fill
    (following ``cellsOnCell``) propagated starting from the seed points and 0
    elsewhere

The command-line tool takes the following arguments:

.. code-block::

    $ compute_mpas_flood_fill_mask --help
    usage: compute_mpas_flood_fill_mask [-h] -m MESH_FILE_NAME -g
                                        GEOJSON_FILE_NAME -o MASK_FILE_NAME

    optional arguments:
      -h, --help            show this help message and exit
      -m MESH_FILE_NAME, --mesh_file_name MESH_FILE_NAME
                            An MPAS mesh file
      -g GEOJSON_FILE_NAME, --geojson_file_name GEOJSON_FILE_NAME
                            An Geojson file containing points at which to start
                            the flood fill
      -o MASK_FILE_NAME, --mask_file_name MASK_FILE_NAME
                            An output MPAS region masks file


Computing Lon/Lat Region Masks
------------------------------

The function :py:func:`mpas_tools.mesh.mask.compute_lon_lat_region_masks()`
or the ``compute_lon_lat_region_masks`` command-line tool compute region masks
on a longitude/latitude grid but are otherwise functionally very similar to
the corresponding tools for compute MPAS region masks. The major difference is
that 1D arrays of longitude and latitude are provided instead of an MPAS mesh
dataset.  There is no argument equivalent to the mask type for MPAS meshes.
Instead, mask values are given at each point on the 2D longitude/latitude grid.
All other arguments serve the same purpose as for the MPAS region mask creation
described above.

The command-line tool takes the following arguments:

.. code-block::

    $ compute_lon_lat_region_masks --help
    usage: compute_lon_lat_region_masks [-h] -i GRID_FILE_NAME [--lon LON]
                                        [--lat LAT] -g GEOJSON_FILE_NAME -o
                                        MASK_FILE_NAME [-c CHUNK_SIZE]
                                        [--show_progress] [-s SUBDIVISION]
                                        [--process_count PROCESS_COUNT]
                                        [--multiprocessing_method MULTIPROCESSING_METHOD]

    optional arguments:
      -h, --help            show this help message and exit
      -i GRID_FILE_NAME, --grid_file_name GRID_FILE_NAME
                            An input lon/lat grid file
      --lon LON             The name of the longitude coordinate
      --lat LAT             The name of the latitude coordinate
      -g GEOJSON_FILE_NAME, --geojson_file_name GEOJSON_FILE_NAME
                            An Geojson file containing mask regions
      -o MASK_FILE_NAME, --mask_file_name MASK_FILE_NAME
                            An output MPAS region masks file
      -c CHUNK_SIZE, --chunk_size CHUNK_SIZE
                            The number of grid points that are processed in one
                            operation
      --show_progress       Whether to show a progress bar
      -s SUBDIVISION, --subdivision SUBDIVISION
                            A threshold in degrees (lon or lat) above which the
                            mask region will be subdivided into smaller polygons
                            for faster intersection checking
      --process_count PROCESS_COUNT
                            The number of processes to use to compute masks. The
                            default is to use all available cores
      --multiprocessing_method MULTIPROCESSING_METHOD
                            The multiprocessing method use for python mask
                            creation ('fork', 'spawn' or 'forkserver')


.. _cull_mpas_dataset:

Culling MPAS Datasets
=====================

The tools described in :ref:`cell_culler` can be used to create a culled
horizontal MPAS mesh.  Once a culled MPAS mesh has been created, an MPAS
dataset on the unculled mesh can be cropped to the culled mesh using the
the :py:func:`mpas_tools.mesh.cull.cull_dataset()` or
:py:func:`mpas_tools.mesh.cull.write_culled_dataset()` functions.  These
functions take a dataset (or filename) to crop as well as datasets (or
filenames) for the unculled and culled horizontal MPAS meshes.  They return
(or write out) the culled version of the data set.  Fields that exist in
the culled horizonal mesh are copied from the culled mesh, rather than cropped
from the dataset.  This because we wish to keep the cropped horizontal mesh
exactly as it was produced by the culling tool, which may not correspond to
a cropped version of the field from the original mesh.  For example, fields
are reindexed during culling and coordinates are recomputed.

It may be useful to compute and store the maps from cells, edges and vertices
on the culled mesh back to the unculled mesh for reuse.  This can be
accomplished by calling the :py:func:`mpas_tools.mesh.cull.map_culled_to_base()`
or :py:func:`mpas_tools.mesh.cull.write_map_culled_to_base()` functions.

An example workflow that culls out ice-shelf cavities from an MPAS-Ocean
initial condition might look like the following.  In this case the file
``culled_mesh.nc`` is a mesh where land (and the grounded portion of the
ice sheet) has been removed but where ice-shelf cavities are still present.
It serves as the "base" mesh for the purposes of this example.
``culled_mesh_no_isc.nc`` is created (if it doesn't already exist) with the
ice-shelf cavities removed as well, so it is the "culled" mesh in this example.
We store the mapping betwen the two horizontal meshes in
``no_isc_to_culled_map.nc`` in case we want to resue it later.  The initial
condition is read from ``initial_state.nc`` and the culled version is written
to ``initial_state_no_isc.nc``:

.. code-block:: python

    import os

    import xarray as xr

    from mpas_tools.io import write_netcdf
    from mpas_tools.mesh.conversion import cull
    from mpas_tools.mesh.cull import write_map_culled_to_base, write_culled_dataset
    from mpas_tools.logging import LoggingContext


    in_filename = 'initial_state.nc'
    out_filename = 'initial_state_no_isc.nc'
    base_mesh_filename = 'culled_mesh.nc'
    culled_mesh_filename = 'culled_mesh_no_isc.nc'
    map_filename = 'no_isc_to_culled_map.nc'

    if not os.path.exists(culled_mesh_filename):
        ds_culled_mesh = xr.open_dataset(base_mesh_filename)
        ds_init = xr.open_dataset(in_filename)
        ds_culled_mesh['cullCell'] = ds_init.landIceMask
        ds_culled_mesh_no_isc = cull(ds_culled_mesh)
        write_netcdf(ds_culled_mesh_no_isc, culled_mesh_filename)

    if not os.path.exists(map_filename):
        write_map_culled_to_base(base_mesh_filename=base_mesh_filename,
                                 culled_mesh_filename=culled_mesh_filename,
                                 out_filename=map_filename)

    with LoggingContext('test') as logger:
        write_culled_dataset(in_filename=in_filename, out_filename=out_filename,
                             base_mesh_filename=base_mesh_filename,
                             culled_mesh_filename=culled_mesh_filename,
                             map_culled_to_base_filename=map_filename,
                             logger=logger)

.. _merge_split:

Merging and Splitting
=====================

In order to support running
`MPAS-Albany Land Ice (MALI) <https://github.com/MPAS-Dev/MPAS-Model/tree/landice/develop>`_
with both Greenland and Antarctica at the same time, tools have been added to
support merging and splitting MPAS meshes.

Merging two meshes can be accomplished with
:py:func:`mpas_tools.merge_grids.merge_grids()`:

.. code-block:: python

    from mpas_tools.translate import translate
    from mpas_tools.merge_grids import merge_grids
    from mpas_tools.planar_hex import make_planar_hex_mesh
    from mpas_tools.io import write_netcdf


    dsMesh1 = make_planar_hex_mesh(nx=10, ny=10, dc=1000., nonperiodic_x=True,
                                   nonperiodic_y=True)

    dsMesh2 = make_planar_hex_mesh(nx=10, ny=10, dc=1000., nonperiodic_x=True,
                                   nonperiodic_y=True)

    translate(dsMesh2, xOffset=20000., yOffset=0.)

    write_netcdf(dsMesh1, 'mesh1.nc')
    write_netcdf(dsMesh2, 'mesh2.nc')

    merge_grids(infile1='mesh1.nc', infile2='mesh2.nc',
                outfile='merged_mesh.nc')

Typically, it will only make sense to merge non-periodic meshes in this way.

Later, perhaps during analysis or visualization, it can be useful to split
apart the merged meshes.  This can be done with
:py:func:`mpas_tools.split_grids.split_grids()`

.. code-block:: python

    from mpas_tools.translate import translate
    from mpas_tools.split_grids import split_grids
    from mpas_tools.planar_hex import make_planar_hex_mesh
    from mpas_tools.io import write_netcdf


    dsMesh1 = make_planar_hex_mesh(nx=10, ny=10, dc=1000., nonperiodic_x=True,
                                   nonperiodic_y=True)

    dsMesh2 = make_planar_hex_mesh(nx=10, ny=10, dc=1000., nonperiodic_x=True,
                                   nonperiodic_y=True)

    translate(dsMesh2, xOffset=20000., yOffset=0.)

    write_netcdf(dsMesh1, 'mesh1.nc')
    write_netcdf(dsMesh2, 'mesh2.nc')


    split_grids(infile='merged_mesh.nc', outfile1='split_mesh1.nc',
                outfile='split_mesh2.nc')

Merging meshes can also be accomplished with the ``merge_grids`` command-line
tool:

.. code-block:: none

    $ merge_grids --help

    usage: merge_grids [-h] [-o FILENAME] FILENAME1 FILENAME2

    Tool to merge 2 MPAS non-contiguous meshes together into a single file

    positional arguments:
      FILENAME1    File name for first mesh to merge
      FILENAME2    File name for second mesh to merge

    optional arguments:
      -h, --help   show this help message and exit
      -o FILENAME  The merged mesh file

Similarly, ``split_grids`` can be used to to split meshes:

.. code-block:: none

    $ split_grids --help

    usage: split_grids [-h] [-1 FILENAME] [-2 FILENAME] [--nCells NCELLS]
                       [--nEdges NEDGES] [--nVertices NVERTICES]
                       [--maxEdges MAXEDGES1 MAXEDGES2]
                       MESHFILE

    Tool to split 2 previously merged MPAS non-contiguous meshes into separate files.
    Typical usage is:
        split_grids.py -1 outfile1.nc -2 outfile2.nc infile
    The optional arguments for nCells, nEdges, nVertices, and maxEdges should
    generally not be required as this information is saved in the combined mesh file
    as global attributes by the merge_grids.py script.

    positional arguments:
      MESHFILE              Mesh file to split

    optional arguments:
      -h, --help            show this help message and exit
      -1 FILENAME, --outfile1 FILENAME
                            File name for first mesh output
                            (default: mesh1.nc)
      -2 FILENAME, --outfile2 FILENAME
                            File name for second mesh output
                            (default: mesh2.nc)
      --nCells NCELLS       The number of cells in the first mesh
                            (default: the value specified in MESHFILE global attribute merge_point)
      --nEdges NEDGES       The number of edges in the first mesh
                            (default: the value specified in MESHFILE global attribute merge_point)
      --nVertices NVERTICES
                            The number of vertices in the first mesh
                            (default: the value specified in MESHFILE global attribute merge_point)
      --maxEdges MAXEDGES1 MAXEDGES2
                            The number of maxEdges in each mesh
                            (default: the value specified in MESHFILE global attribute merge_point
                                  OR: will use MESHFILE maxEdges dimension and assume same for both)


.. _mesh_translation:

Translation
===========

A planar mesh can be translated in x, y or both by calling
:py:func:`mpas_tools.translate.translate()`:

.. code-block:: python

    from mpas_tools.translate import translate
    from mpas_tools.planar_hex import make_planar_hex_mesh

    dsMesh = make_planar_hex_mesh(nx=10, ny=20, dc=1000., nonperiodic_x=False,
                                  nonperiodic_y=False)

    translate(dsMesh, xOffset=1000., yOffset=2000.)

This creates a periodic, planar mesh and then translates it by 1 km in x and
2 km in y.

.. note::

    All the functions in the ``mpas_tools.translate`` module modify the mesh
    inplace, rather than returning a new ``xarray.Dataset`` object.  This is
    in contrast to typical ``xarray`` functions and methods.


A mesh can be translated so that its center is at ``x = 0.``, ``y = 0.`` with
the function :py:func:`mpas_tools.translate.center()`:

.. code-block:: python

    from mpas_tools.translate import center
    from mpas_tools.planar_hex import make_planar_hex_mesh

    dsMesh = make_planar_hex_mesh(nx=10, ny=20, dc=1000., nonperiodic_x=False,
                                  nonperiodic_y=False)

    center(dsMesh)

A mesh can be translated so its center matches the center of another mesh by
using :py:func:`mpas_tools.translate.center_on_mesh()`:

.. code-block:: python

    from mpas_tools.translate import center_on_mesh
    from mpas_tools.planar_hex import make_planar_hex_mesh

    dsMesh1 = make_planar_hex_mesh(nx=10, ny=20, dc=1000., nonperiodic_x=False,
                                   nonperiodic_y=False)

    dsMesh2 = make_planar_hex_mesh(nx=20, ny=40, dc=2000., nonperiodic_x=False,
                                   nonperiodic_y=False)

    center_on_mesh(dsMesh2, dsMesh1)

In this example, the coordinates of ``dsMesh2`` are altered so its center
matches that of ``dsMesh1``.

The functionality of all three of these functions is also available via the
``translate_planar_grid`` command-line tool:

.. code-block:: none

    $ translate_planar_grid --help

    == Gathering information.  (Invoke with --help for more details. All arguments are optional)
    Usage: translate_planar_grid [options]

    This script translates the coordinate system of the planar MPAS mesh specified
    with the -f flag.  There are 3 possible methods to choose from: 1) shift the
    origin to the center of the domain 2) arbirary shift in x and/or y 3) shift to
    the center of the domain described in a separate file

    Options:
      -h, --help            show this help message and exit
      -f FILENAME, --file=FILENAME
                            MPAS planar grid file name. [default: grid.nc]
      -d FILENAME, --datafile=FILENAME
                            data file name to which to match the domain center of.
                            Uses xCell,yCell or, if those fields do not exist,
                            will secondly try x1,y1 fields.
      -x SHIFT_VALUE        user-specified shift in the x-direction. [default:
                            0.0]
      -y SHIFT_VALUE        user-specified shift in the y-direction. [default:
                            0.0]
      -c                    shift so origin is at center of domain [default:
                            False]


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

A command-line tools is also available for this purpose:

.. code-block:: none

    $ scrip_from_mpas --help
    == Gathering information.  (Invoke with --help for more details. All arguments are optional)
    Usage: scrip_from_mpas [options]

    This script takes an MPAS grid file and generates a SCRIP grid file.

    Options:
      -h, --help            show this help message and exit
      -m FILENAME, --mpas=FILENAME
                            MPAS grid file name used as input. [default: grid.nc]
      -s FILENAME, --scrip=FILENAME
                            SCRIP grid file to output. [default: scrip.nc]
      -l, --landice         If flag is on, landice masks will be computed and
                            used.
