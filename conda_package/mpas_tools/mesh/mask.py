import argparse
from functools import partial

import numpy
import progressbar
import shapely.geometry
import shapely.ops
import xarray as xr
from geometric_features import read_feature_collection
from igraph import Graph
from scipy.spatial import KDTree
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, box
from shapely.strtree import STRtree

from mpas_tools.cime.constants import constants
from mpas_tools.io import default_nchar, write_netcdf
from mpas_tools.logging import LoggingContext
from mpas_tools.parallel import create_pool
from mpas_tools.transects import (
    cartesian_to_lon_lat,
    lon_lat_to_cartesian,
    subdivide_great_circle,
)


def compute_mpas_region_masks(
    dsMesh,
    fcMask,
    maskTypes=('cell', 'vertex'),
    logger=None,
    pool=None,
    chunkSize=1000,
    showProgress=False,
    subdivisionThreshold=30.0,
):
    """
    Use shapely and processes to create a set of masks from a feature
    collection made up of regions (polygons)

    Parameters
    ----------
    dsMesh : xarray.Dataset
        An MPAS mesh on which the masks should be created

    fcMask : geometric_features.FeatureCollection
        A feature collection containing features to use to create the mask

    maskTypes : tuple of {'cell', 'edge', 'vertex'}, optional
        Which type(s) of masks to make.  Masks are created based on whether
        the latitude and longitude associated with each of these locations
        (e.g. ``dsMesh.latCell`` and ``dsMesh.lonCell`` for ``'cell'``) are
        inside or outside of the regions in ``fcMask``.

    logger : logging.Logger, optional
        A logger for the output if not stdout

    pool : multiprocessing.Pool, optional
        A pool for performing multiprocessing

    chunkSize : int, optional
        The number of cells, vertices or edges that are processed in one
        operation.  Experimentation has shown that 1000 is a reasonable
        compromise between dividing the work into sufficient subtasks to
        distribute the load and having sufficient work for each thread.

    showProgress : bool, optional
        Whether to show a progress bar

    subdivisionThreshold : float, optional
        A threshold in degrees (lon or lat) above which the mask region will
        be subdivided into smaller polygons for faster intersection checking

    Returns
    -------
    dsMask : xarray.Dataset
        The masks

    """

    suffixes = {'cell': 'Cell', 'edge': 'Edge', 'vertex': 'Vertex'}
    dims = {'cell': 'nCells', 'edge': 'nEdges', 'vertex': 'nVertices'}

    dsMasks = xr.Dataset()

    for maskType in maskTypes:
        suffix = suffixes[maskType]
        dim = dims[maskType]
        lonName = f'lon{suffix}'
        latName = f'lat{suffix}'
        lat = numpy.rad2deg(dsMesh[latName].values)

        # transform longitudes to [-180, 180)
        lon = (
            numpy.mod(numpy.rad2deg(dsMesh[lonName].values) + 180.0, 360.0)
            - 180.0
        )

        if logger is not None:
            logger.info(f'  Computing {maskType} masks:')

        # create shapely geometry for lon and lat
        points = [shapely.geometry.Point(x, y) for x, y in zip(lon, lat)]
        regionNames, masks, properties = _compute_region_masks(
            fcMask,
            points,
            logger,
            pool,
            chunkSize,
            showProgress,
            subdivisionThreshold,
        )

        nPoints = len(points)

        if logger is not None:
            logger.info('  Adding masks to dataset...')
        nRegions = len(regionNames)
        # create a new data array for masks
        masksVarName = f'region{suffix}Masks'
        dsMasks[masksVarName] = (
            (dim, 'nRegions'),
            numpy.zeros((nPoints, nRegions), dtype=numpy.int32),
        )

        for index in range(nRegions):
            mask = masks[index]
            dsMasks[masksVarName][:, index] = numpy.array(
                mask, dtype=numpy.int32
            )

        properties['regionNames'] = regionNames
        _add_properties(
            ds=dsMasks,
            properties=properties,
            dim='nRegions',
        )

    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def entry_point_compute_mpas_region_masks():
    """Entry point for ``compute_mpas_region_masks()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m',
        '--mesh_file_name',
        dest='mesh_file_name',
        type=str,
        required=True,
        help='An MPAS mesh file',
    )
    parser.add_argument(
        '-g',
        '--geojson_file_name',
        dest='geojson_file_name',
        type=str,
        required=True,
        help='An Geojson file containing mask regions',
    )
    parser.add_argument(
        '-o',
        '--mask_file_name',
        dest='mask_file_name',
        type=str,
        required=True,
        help='An output MPAS region masks file',
    )
    parser.add_argument(
        '-t',
        '--mask_types',
        nargs='+',
        dest='mask_types',
        type=str,
        help='Which type(s) of masks to make: cell, edge or '
        'vertex.  Default is cell and vertex.',
    )
    parser.add_argument(
        '-c',
        '--chunk_size',
        dest='chunk_size',
        type=int,
        default=1000,
        help='The number of cells, vertices or edges that are '
        'processed in one operation',
    )
    parser.add_argument(
        '--show_progress',
        dest='show_progress',
        action='store_true',
        help='Whether to show a progress bar',
    )
    parser.add_argument(
        '-s',
        '--subdivision',
        dest='subdivision',
        type=float,
        default=30.0,
        help='A threshold in degrees (lon or lat) above which '
        'the mask region will be subdivided into smaller '
        'polygons for faster intersection checking',
    )
    parser.add_argument(
        '--process_count',
        required=False,
        dest='process_count',
        type=int,
        help='The number of processes to use to compute masks.  The '
        'default is to use all available cores',
    )
    parser.add_argument(
        '--multiprocessing_method',
        dest='multiprocessing_method',
        default='forkserver',
        help='The multiprocessing method use for python mask creation '
        "('fork', 'spawn' or 'forkserver')",
    )
    parser.add_argument(
        '--format', dest='format', type=str, help='NetCDF file format'
    )
    parser.add_argument(
        '--engine', dest='engine', type=str, help='NetCDF output engine'
    )
    args = parser.parse_args()

    dsMesh = xr.open_dataset(
        args.mesh_file_name, decode_cf=False, decode_times=False
    )
    fcMask = read_feature_collection(args.geojson_file_name)

    pool = create_pool(
        process_count=args.process_count, method=args.multiprocessing_method
    )

    if args.mask_types is None:
        args.mask_types = ('cell', 'vertex')

    with LoggingContext('compute_mpas_region_masks') as logger:
        dsMasks = compute_mpas_region_masks(
            dsMesh=dsMesh,
            fcMask=fcMask,
            maskTypes=args.mask_types,
            logger=logger,
            pool=pool,
            chunkSize=args.chunk_size,
            showProgress=args.show_progress,
            subdivisionThreshold=args.subdivision,
        )

        write_netcdf(
            dsMasks,
            args.mask_file_name,
            format=args.format,
            engine=args.engine,
            logger=logger,
        )


def compute_mpas_transect_masks(
    dsMesh,
    fcMask,
    earthRadius,
    maskTypes=('cell', 'edge', 'vertex'),
    logger=None,
    pool=None,
    chunkSize=1000,
    showProgress=False,
    subdivisionResolution=10e3,
    addEdgeSign=False,
):
    """
    Use shapely and processes to create a set of masks from a feature
    collection made up of transects (line strings)

    Parameters
    ----------
    dsMesh : xarray.Dataset
        An MPAS mesh on which the masks should be created

    fcMask : geometric_features.FeatureCollection
        A feature collection containing features to use to create the mask

    earthRadius : float
        The radius of the earth in meters

    maskTypes : tuple of {'cell', 'edge', 'vertex'}, optional
        Which type(s) of masks to make.  Masks are created based on whether
        the latitude and longitude associated with each of these locations
        (e.g. ``dsMesh.latCell`` and ``dsMesh.lonCell`` for ``'cell'``) are
        inside or outside of the transects in ``fcMask``.

    logger : logging.Logger, optional
        A logger for the output if not stdout

    pool : multiprocessing.Pool, optional
        A pool for performing multiprocessing

    chunkSize : int, optional
        The number of cells, vertices or edges that are processed in one
        operation.  Experimentation has shown that 1000 is a reasonable
        compromise between dividing the work into sufficient subtasks to
        distribute the load and having sufficient work for each thread.

    showProgress : bool, optional
        Whether to show a progress bar

    subdivisionResolution : float, optional
        The maximum resolution (in meters) of segments in a transect.  If a
        transect is too coarse, it will be subdivided.  Pass ``None`` for no
        subdivision.

    addEdgeSign : bool, optional
        Whether to add the ``edgeSign`` variable, which requires significant
        extra computation

    Returns
    -------
    dsMask : xarray.Dataset
        The masks

    """

    suffixes = {'cell': 'Cell', 'edge': 'Edge', 'vertex': 'Vertex'}
    dims = {'cell': 'nCells', 'edge': 'nEdges', 'vertex': 'nVertices'}

    dsMasks = xr.Dataset()

    for maskType in maskTypes:
        suffix = suffixes[maskType]
        dim = dims[maskType]

        if logger is not None:
            logger.info(f'  Computing {maskType} masks:')

        polygons, nPolygons, duplicatePolygons = _get_polygons(
            dsMesh, maskType
        )
        transectNames, masks, properties, shapes = _compute_transect_masks(
            fcMask,
            polygons,
            logger,
            pool,
            chunkSize,
            showProgress,
            subdivisionResolution,
            earthRadius,
        )

        if logger is not None:
            if addEdgeSign and maskType == 'edge':
                logger.info('  Adding masks and edge signs to dataset...')
            else:
                logger.info('  Adding masks to dataset...')
        nTransects = len(transectNames)
        # create a new data array for masks
        masksVarName = f'transect{suffix}Masks'
        dsMasks[masksVarName] = (
            (dim, 'nTransects'),
            numpy.zeros((nPolygons, nTransects), dtype=numpy.int32),
        )

        if addEdgeSign and maskType == 'edge':
            dsMasks['transectEdgeMaskSigns'] = (
                (dim, 'nTransects'),
                numpy.zeros((nPolygons, nTransects), dtype=numpy.int32),
            )

        for index in range(nTransects):
            maskAndDuplicates = masks[index]
            mask = maskAndDuplicates[0:nPolygons]

            mask[duplicatePolygons] = numpy.logical_or(
                mask[duplicatePolygons], maskAndDuplicates[nPolygons:]
            )
            dsMasks[masksVarName][:, index] = numpy.array(
                mask, dtype=numpy.int32
            )

            if addEdgeSign and maskType == 'edge':
                print(transectNames[index])
                dsMasks['transectEdgeMaskSigns'][:, index] = (
                    _compute_edge_sign(dsMesh, mask, shapes[index])
                )

        properties['transectNames'] = transectNames
        _add_properties(
            ds=dsMasks,
            properties=properties,
            dim='nTransects',
        )

    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def entry_point_compute_mpas_transect_masks():
    """Entry point for ``compute_mpas_transect_masks()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m',
        '--mesh_file_name',
        dest='mesh_file_name',
        type=str,
        required=True,
        help='An MPAS mesh file',
    )
    parser.add_argument(
        '-g',
        '--geojson_file_name',
        dest='geojson_file_name',
        type=str,
        required=True,
        help='An Geojson file containing transects',
    )
    parser.add_argument(
        '-o',
        '--mask_file_name',
        dest='mask_file_name',
        type=str,
        required=True,
        help='An output MPAS transect masks file',
    )
    parser.add_argument(
        '-t',
        '--mask_types',
        nargs='+',
        dest='mask_types',
        type=str,
        help='Which type(s) of masks to make: cell, edge or '
        'vertex.  Default is cell, edge and vertex.',
    )
    parser.add_argument(
        '-c',
        '--chunk_size',
        dest='chunk_size',
        type=int,
        default=1000,
        help='The number of cells, vertices or edges that are '
        'processed in one operation',
    )
    parser.add_argument(
        '--show_progress',
        dest='show_progress',
        action='store_true',
        help='Whether to show a progress bar',
    )
    parser.add_argument(
        '-s',
        '--subdivision',
        dest='subdivision',
        type=float,
        help='The maximum resolution (in meters) of segments '
        'in a transect.  If a transect is too coarse, it '
        'will be subdivided.  Default is no subdivision.',
    )
    parser.add_argument(
        '--process_count',
        required=False,
        dest='process_count',
        type=int,
        help='The number of processes to use to compute masks.  The '
        'default is to use all available cores',
    )
    parser.add_argument(
        '--multiprocessing_method',
        dest='multiprocessing_method',
        default='forkserver',
        help='The multiprocessing method use for python mask creation '
        "('fork', 'spawn' or 'forkserver')",
    )
    parser.add_argument(
        '--add_edge_sign',
        dest='add_edge_sign',
        action='store_true',
        help='Whether to add the transectEdgeMaskSigns variable',
    )
    parser.add_argument(
        '--format', dest='format', type=str, help='NetCDF file format'
    )
    parser.add_argument(
        '--engine', dest='engine', type=str, help='NetCDF output engine'
    )
    args = parser.parse_args()

    dsMesh = xr.open_dataset(
        args.mesh_file_name, decode_cf=False, decode_times=False
    )
    fcMask = read_feature_collection(args.geojson_file_name)

    pool = create_pool(
        process_count=args.process_count, method=args.multiprocessing_method
    )

    if args.mask_types is None:
        args.mask_types = ('cell', 'edge', 'vertex')

    earth_radius = constants['SHR_CONST_REARTH']

    with LoggingContext('compute_mpas_transect_masks') as logger:
        dsMasks = compute_mpas_transect_masks(
            dsMesh=dsMesh,
            fcMask=fcMask,
            earthRadius=earth_radius,
            maskTypes=args.mask_types,
            logger=logger,
            pool=pool,
            chunkSize=args.chunk_size,
            showProgress=args.show_progress,
            subdivisionResolution=args.subdivision,
            addEdgeSign=args.add_edge_sign,
        )

        write_netcdf(
            dsMasks,
            args.mask_file_name,
            format=args.format,
            engine=args.engine,
            logger=logger,
        )


def compute_mpas_flood_fill_mask(
    dsMesh, fcSeed, daGrow=None, logger=None, workers=-1
):
    """
    Flood fill from the given set of seed points to create a contiguous mask.
    The flood fill operates using cellsOnCell, starting from the cells
    whose centers are closest to the seed points.

    Parameters
    ----------
    dsMesh : xarray.Dataset
        An MPAS mesh on which the masks should be created

    fcSeed : geometric_features.FeatureCollection
        A feature collection containing points at which to start the flood fill

    daGrow : xarray.DataArray, optional
        A data array of size ``nCells`` with a mask that is 1 anywhere the
        flood fill is allowed to grow.  The default is that the mask is all
        ones.

    logger : logging.Logger, optional
        A logger for the output if not stdout

    workers : int, optional
        The number of threads used for finding nearest neighbors.  The default
        is all available threads (``workers=-1``)

    Returns
    -------
    dsMask : xarray.Dataset
        The masks

    """

    dsMasks = xr.Dataset()

    lat = numpy.rad2deg(dsMesh.latCell.values)

    # transform longitudes to [-180, 180)
    lon = (
        numpy.mod(numpy.rad2deg(dsMesh.lonCell.values) + 180.0, 360.0) - 180.0
    )

    if logger is not None:
        logger.info('  Computing flood fill mask on cells:')

    seedMask = _compute_seed_mask(fcSeed, lon, lat, workers)

    cellsOnCell = dsMesh.cellsOnCell.values - 1

    if daGrow is not None:
        growMask = daGrow.values
    else:
        growMask = numpy.ones(dsMesh.sizes['nCells'])

    seedMask = _flood_fill_mask(seedMask, growMask, cellsOnCell)

    if logger is not None:
        logger.info('  Adding masks to dataset...')
    # create a new data array for the mask
    masksVarName = 'cellSeedMask'
    dsMasks[masksVarName] = (
        ('nCells',),
        numpy.array(seedMask, dtype=numpy.int32),
    )

    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def entry_point_compute_mpas_flood_fill_mask():
    """Entry point for ``compute_mpas_flood_fill_mask()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m',
        '--mesh_file_name',
        dest='mesh_file_name',
        type=str,
        required=True,
        help='An MPAS mesh file',
    )
    parser.add_argument(
        '-g',
        '--geojson_file_name',
        dest='geojson_file_name',
        type=str,
        required=True,
        help='An Geojson file containing points at which to '
        'start the flood fill',
    )
    parser.add_argument(
        '-o',
        '--mask_file_name',
        dest='mask_file_name',
        type=str,
        required=True,
        help='An output MPAS region masks file',
    )
    parser.add_argument(
        '--format', dest='format', type=str, help='NetCDF file format'
    )
    parser.add_argument(
        '--engine', dest='engine', type=str, help='NetCDF output engine'
    )
    args = parser.parse_args()

    dsMesh = xr.open_dataset(
        args.mesh_file_name, decode_cf=False, decode_times=False
    )
    fcSeed = read_feature_collection(args.geojson_file_name)

    with LoggingContext('compute_mpas_flood_fill_mask') as logger:
        dsMasks = compute_mpas_flood_fill_mask(
            dsMesh=dsMesh, fcSeed=fcSeed, logger=logger
        )

        write_netcdf(
            dsMasks,
            args.mask_file_name,
            format=args.format,
            engine=args.engine,
            logger=logger,
        )


def compute_lon_lat_region_masks(
    lon,
    lat,
    fcMask,
    logger=None,
    pool=None,
    chunkSize=1000,
    showProgress=False,
    subdivisionThreshold=30.0,
):
    """
    Use shapely and processes to create a set of masks from a feature
    collection made up of regions (polygons) on a tensor lon/lat grid

    Parameters
    ----------
    lon : numpy.ndarray
        A 1D array of longitudes in degrees between -180 and 180

    lat : numpy.ndarray
        A 1D array of latitudes in degrees between -90 and 90

    fcMask : geometric_features.FeatureCollection
        A feature collection containing features to use to create the mask

    logger : logging.Logger, optional
        A logger for the output if not stdout

    pool : multiprocessing.Pool, optional
        A pool for performing multiprocessing

    chunkSize : int, optional
        The number of cells, vertices or edges that are processed in one
        operation.  Experimentation has shown that 1000 is a reasonable
        compromise between dividing the work into sufficient subtasks to
        distribute the load and having sufficient work for each thread.

    showProgress : bool, optional
        Whether to show a progress bar

    subdivisionThreshold : float, optional
        A threshold in degrees (lon or lat) above which the mask region will
        be subdivided into smaller polygons for faster intersection checking

    Returns
    -------
    dsMask : xarray.Dataset
        The masks

    """

    dsMasks = xr.Dataset()

    # make sure lon is between -180 and 180
    lon = numpy.mod(lon + 180.0, 360.0) - 180.0

    Lon, Lat = numpy.meshgrid(lon, lat)

    shape = Lon.shape

    Lon = Lon.ravel()
    Lat = Lat.ravel()

    # create shapely geometry for lon and lat
    points = [shapely.geometry.Point(x, y) for x, y in zip(Lon, Lat)]
    regionNames, masks, properties = _compute_region_masks(
        fcMask,
        points,
        logger,
        pool,
        chunkSize,
        showProgress,
        subdivisionThreshold,
    )

    nlon = len(lon)
    nlat = len(lat)

    if logger is not None:
        logger.info('  Adding masks to dataset...')
    nRegions = len(regionNames)
    # create a new data array for masks
    masksVarName = 'regionMasks'
    dsMasks[masksVarName] = (
        ('lat', 'lon', 'nRegions'),
        numpy.zeros((nlat, nlon, nRegions), dtype=numpy.int32),
    )

    for index in range(nRegions):
        mask = masks[index]
        dsMasks[masksVarName][:, :, index] = numpy.array(
            mask.reshape(shape), dtype=numpy.int32
        )

    properties['regionNames'] = regionNames
    _add_properties(
        ds=dsMasks,
        properties=properties,
        dim='nRegions',
    )

    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def entry_point_compute_lon_lat_region_masks():
    """Entry point for ``compute_lon_lat_region_masks()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--grid_file_name',
        dest='grid_file_name',
        type=str,
        required=True,
        help='An input lon/lat grid file',
    )
    parser.add_argument(
        '--lon',
        dest='lon',
        default='lon',
        type=str,
        help='The name of the longitude coordinate',
    )
    parser.add_argument(
        '--lat',
        dest='lat',
        default='lat',
        type=str,
        help='The name of the latitude coordinate',
    )
    parser.add_argument(
        '-g',
        '--geojson_file_name',
        dest='geojson_file_name',
        type=str,
        required=True,
        help='An Geojson file containing mask regions',
    )
    parser.add_argument(
        '-o',
        '--mask_file_name',
        dest='mask_file_name',
        type=str,
        required=True,
        help='An output MPAS region masks file',
    )
    parser.add_argument(
        '-c',
        '--chunk_size',
        dest='chunk_size',
        type=int,
        default=1000,
        help='The number of grid points that are processed in one operation',
    )
    parser.add_argument(
        '--show_progress',
        dest='show_progress',
        action='store_true',
        help='Whether to show a progress bar',
    )
    parser.add_argument(
        '-s',
        '--subdivision',
        dest='subdivision',
        type=float,
        default=30.0,
        help='A threshold in degrees (lon or lat) above which '
        'the mask region will be subdivided into smaller '
        'polygons for faster intersection checking',
    )
    parser.add_argument(
        '--process_count',
        required=False,
        dest='process_count',
        type=int,
        help='The number of processes to use to compute masks.  The '
        'default is to use all available cores',
    )
    parser.add_argument(
        '--multiprocessing_method',
        dest='multiprocessing_method',
        default='forkserver',
        help='The multiprocessing method use for python mask creation '
        "('fork', 'spawn' or 'forkserver')",
    )
    parser.add_argument(
        '--format', dest='format', type=str, help='NetCDF file format'
    )
    parser.add_argument(
        '--engine', dest='engine', type=str, help='NetCDF output engine'
    )
    args = parser.parse_args()

    dsGrid = xr.open_dataset(
        args.grid_file_name, decode_cf=False, decode_times=False
    )
    lon = dsGrid[args.lon].values
    lat = dsGrid[args.lat].values

    fcMask = read_feature_collection(args.geojson_file_name)

    pool = create_pool(
        process_count=args.process_count, method=args.multiprocessing_method
    )

    with LoggingContext('compute_lon_lat_region_masks') as logger:
        dsMasks = compute_lon_lat_region_masks(
            lon=lon,
            lat=lat,
            fcMask=fcMask,
            logger=logger,
            pool=pool,
            chunkSize=args.chunk_size,
            showProgress=args.show_progress,
            subdivisionThreshold=args.subdivision,
        )

    write_netcdf(
        dsMasks,
        args.mask_file_name,
        format=args.format,
        engine=args.engine,
        logger=logger,
    )


def compute_projection_grid_region_masks(
    lon,
    lat,
    fcMask,
    logger=None,
    pool=None,
    chunkSize=1000,
    showProgress=False,
    subdivisionThreshold=30.0,
    xdim='x',
    ydim='y',
):
    """
    Use shapely and processes to create a set of masks from a feature
    collection made up of regions (polygons) on a projection grid such as
    a polar-stereographic grid.

    Parameters
    ----------
    lon : numpy.ndarray
        A 2D array of longitudes in degrees between -180 and 180

    lat : numpy.ndarray
        A 2D array of latitudes in degrees between -90 and 90

    fcMask : geometric_features.FeatureCollection
        A feature collection containing features to use to create the mask

    logger : logging.Logger, optional
        A logger for the output if not stdout

    pool : multiprocessing.Pool, optional
        A pool for performing multiprocessing

    chunkSize : int, optional
        The number of cells, vertices or edges that are processed in one
        operation.  Experimentation has shown that 1000 is a reasonable
        compromise between dividing the work into sufficient subtasks to
        distribute the load and having sufficient work for each thread.

    showProgress : bool, optional
        Whether to show a progress bar

    subdivisionThreshold : float, optional
        A threshold in degrees (lon or lat) above which the mask region will
        be subdivided into smaller polygons for faster intersection checking

    xdim : str, optional
        The name of the x dimension

    ydim : str, optional
        The name of the y dimension

    Returns
    -------
    dsMask : xarray.Dataset
        The masks

    """

    dsMasks = xr.Dataset()

    # make sure -180 <= lon < 180
    lon = numpy.mod(lon + 180.0, 360.0) - 180.0

    ny, nx = lon.shape

    # create shapely geometry for lon and lat
    points = [
        shapely.geometry.Point(x, y) for x, y in zip(lon.ravel(), lat.ravel())
    ]
    regionNames, masks, properties = _compute_region_masks(
        fcMask,
        points,
        logger,
        pool,
        chunkSize,
        showProgress,
        subdivisionThreshold,
    )

    if logger is not None:
        logger.info('  Adding masks to dataset...')
    nRegions = len(regionNames)
    # create a new data array for masks
    masksVarName = 'regionMasks'
    dsMasks[masksVarName] = (
        (ydim, xdim, 'nRegions'),
        numpy.zeros((ny, nx, nRegions), dtype=numpy.int32),
    )

    for index in range(nRegions):
        mask = masks[index]
        dsMasks[masksVarName][:, :, index] = numpy.array(
            mask.reshape((ny, nx)), dtype=numpy.int32
        )

    properties['regionNames'] = regionNames
    _add_properties(
        ds=dsMasks,
        properties=properties,
        dim='nRegions',
    )

    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def entry_point_compute_projection_grid_region_masks():
    """Entry point for ``compute_projection_grid_region_masks()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--grid_file_name',
        dest='grid_file_name',
        type=str,
        required=True,
        help='An input lon/lat grid file',
    )
    parser.add_argument(
        '--lon',
        dest='lon',
        default='lon',
        type=str,
        help='The name of the 2D longitude coordinate',
    )
    parser.add_argument(
        '--lat',
        dest='lat',
        default='lat',
        type=str,
        help='The name of the 2D latitude coordinate',
    )
    parser.add_argument(
        '-g',
        '--geojson_file_name',
        dest='geojson_file_name',
        type=str,
        required=True,
        help='An Geojson file containing mask regions',
    )
    parser.add_argument(
        '-o',
        '--mask_file_name',
        dest='mask_file_name',
        type=str,
        required=True,
        help='An output MPAS region masks file',
    )
    parser.add_argument(
        '-c',
        '--chunk_size',
        dest='chunk_size',
        type=int,
        default=1000,
        help='The number of grid points that are processed in one operation',
    )
    parser.add_argument(
        '--show_progress',
        dest='show_progress',
        action='store_true',
        help='Whether to show a progress bar',
    )
    parser.add_argument(
        '-s',
        '--subdivision',
        dest='subdivision',
        type=float,
        default=30.0,
        help='A threshold in degrees (lon or lat) above which '
        'the mask region will be subdivided into smaller '
        'polygons for faster intersection checking',
    )
    parser.add_argument(
        '--process_count',
        required=False,
        dest='process_count',
        type=int,
        help='The number of processes to use to compute masks.  The '
        'default is to use all available cores',
    )
    parser.add_argument(
        '--multiprocessing_method',
        dest='multiprocessing_method',
        default='forkserver',
        help='The multiprocessing method use for python mask creation '
        "('fork', 'spawn' or 'forkserver')",
    )
    parser.add_argument(
        '--format', dest='format', type=str, help='NetCDF file format'
    )
    parser.add_argument(
        '--engine', dest='engine', type=str, help='NetCDF output engine'
    )
    args = parser.parse_args()

    dsGrid = xr.open_dataset(
        args.grid_file_name, decode_cf=False, decode_times=False
    )
    lon = dsGrid[args.lon]
    lat = dsGrid[args.lat]

    ydim, xdim = lon.dims

    fcMask = read_feature_collection(args.geojson_file_name)

    pool = create_pool(
        process_count=args.process_count, method=args.multiprocessing_method
    )

    with LoggingContext('compute_lon_lat_region_masks') as logger:
        dsMasks = compute_projection_grid_region_masks(
            lon=lon.values,
            lat=lat.values,
            fcMask=fcMask,
            logger=logger,
            pool=pool,
            chunkSize=args.chunk_size,
            showProgress=args.show_progress,
            subdivisionThreshold=args.subdivision,
            xdim=xdim,
            ydim=ydim,
        )

    write_netcdf(
        dsMasks,
        args.mask_file_name,
        format=args.format,
        engine=args.engine,
        logger=logger,
    )


def _compute_mask_from_shapes(
    shapes1, shapes2, func, pool, chunkSize, showProgress
):
    """
    If multiprocessing, break shapes2 into chunks and use multiprocessing to
    apply the given function one chunk at a time
    """
    nShapes2 = len(shapes2)
    if pool is None:
        mask = func(shapes1, shapes2)
    else:
        nChunks = int(numpy.ceil(nShapes2 / chunkSize))
        chunks = []
        indices = [0]
        for iChunk in range(nChunks):
            start = iChunk * chunkSize
            end = min((iChunk + 1) * chunkSize, nShapes2)
            chunks.append(shapes2[start:end])
            indices.append(end)

        partial_func = partial(func, shapes1)
        if showProgress:
            widgets = [
                '    ',
                progressbar.Percentage(),
                ' ',
                progressbar.Bar(),
                ' ',
                progressbar.ETA(),
            ]
            bar = progressbar.ProgressBar(
                widgets=widgets, maxval=nChunks
            ).start()
        else:
            bar = None

        mask = numpy.zeros((nShapes2,), bool)
        for iChunk, maskChunk in enumerate(pool.imap(partial_func, chunks)):
            mask[indices[iChunk] : indices[iChunk + 1]] = maskChunk
            if showProgress:
                bar.update(iChunk + 1)
        if showProgress:
            bar.finish()
    return mask


def _add_properties(ds, properties, dim):
    """
    Add properties to the dataset from a dictionary of properties
    """
    nchar = default_nchar
    for name, prop_list in properties.items():
        if name not in ds:
            if isinstance(prop_list[0], str):
                ds[name] = (
                    (dim,),
                    numpy.zeros((len(prop_list),), dtype=f'|S{nchar}'),
                )

                for index, value in enumerate(prop_list):
                    ds[name][index] = value
            else:
                ds[name] = ((dim,), prop_list)


def _get_region_names_and_properties(fc):
    regionNames = []
    for feature in fc.features:
        name = feature['properties']['name']
        regionNames.append(name)

    propertyNames = set()
    for feature in fc.features:
        for propertyName in feature['properties']:
            if propertyName not in [
                'name',
                'author',
                'tags',
                'component',
                'object',
            ]:
                propertyNames.add(propertyName)

    properties = {}
    for propertyName in propertyNames:
        properties[propertyName] = []
        for feature in fc.features:
            if propertyName in feature['properties']:
                propertyVal = feature['properties'][propertyName]
                properties[propertyName].append(propertyVal)
            else:
                properties[propertyName].append('')

    return regionNames, properties


def _compute_region_masks(
    fcMask, points, logger, pool, chunkSize, showProgress, threshold
):
    """
    Build a region mask file from the given mesh and geojson file defining
    a set of regions.
    """

    regionNames, properties = _get_region_names_and_properties(fcMask)

    masks = []

    for feature in fcMask.features:
        name = feature['properties']['name']

        if logger is not None:
            logger.info(f'    {name}')

        shape = shapely.geometry.shape(feature['geometry'])
        shapes = _katana(shape, threshold=threshold)

        mask = _compute_mask_from_shapes(
            shapes1=shapes,
            shapes2=points,
            func=_contains,
            pool=pool,
            chunkSize=chunkSize,
            showProgress=showProgress,
        )

        masks.append(mask)

    return regionNames, masks, properties


def _contains(shapes, points):
    tree = STRtree(points)
    mask = numpy.zeros(len(points), dtype=bool)
    for shape in shapes:
        indicesInShape = tree.query(shape, predicate='covers')
        mask[indicesInShape] = True
    return mask


def _katana(geometry, threshold, count=0, maxcount=250):
    """
    From https://snorfalorpagus.net/blog/2016/03/13/splitting-large-polygons-for-faster-intersections/

    Split a Polygon into two parts across it's shortest dimension

    Copyright (c) 2016, Joshua Arnott

    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    """
    bounds = geometry.bounds
    width = bounds[2] - bounds[0]
    height = bounds[3] - bounds[1]
    if max(width, height) <= threshold or count == maxcount:
        # either the polygon is smaller than the threshold, or the maximum
        # number of recursions has been reached
        return [geometry]
    if height >= width:
        # split left to right
        a = box(bounds[0], bounds[1], bounds[2], bounds[1] + height / 2)
        b = box(bounds[0], bounds[1] + height / 2, bounds[2], bounds[3])
    else:
        # split top to bottom
        a = box(bounds[0], bounds[1], bounds[0] + width / 2, bounds[3])
        b = box(bounds[0] + width / 2, bounds[1], bounds[2], bounds[3])
    result = []
    for d in (
        a,
        b,
    ):
        c = geometry.intersection(d)
        if isinstance(c, GeometryCollection):
            c = c.geoms
        else:
            c = [c]
        for e in c:
            if isinstance(e, (Polygon, MultiPolygon)):
                result.extend(_katana(e, threshold, count + 1, maxcount))
    if count > 0:
        return result
    # convert multipart into singlepart
    final_result = []
    for g in result:
        if isinstance(g, MultiPolygon):
            final_result.extend(g.geoms)
        else:
            final_result.append(g)
    return final_result


def _compute_transect_masks(
    fcMask,
    polygons,
    logger,
    pool,
    chunkSize,
    showProgress,
    subdivisionResolution,
    earthRadius,
):
    """
    Build a transect mask file from the given mesh and geojson file defining
    a set of transects.
    """

    transectNames, properties = _get_region_names_and_properties(fcMask)

    masks = []
    shapes = []

    for feature in fcMask.features:
        name = feature['properties']['name']

        if logger is not None:
            logger.info(f'    {name}')

        geom_type = feature['geometry']['type']
        if geom_type == 'LineString':
            coordinates = [feature['geometry']['coordinates']]
        elif geom_type == 'MultiLineString':
            coordinates = feature['geometry']['coordinates']
        else:
            raise ValueError(f'Unexpected geometry type {geom_type}')

        new_coords = []
        for coords in coordinates:
            if subdivisionResolution is None:
                new_coords.append(coords)
            else:
                lon, lat = zip(*coords)
                x, y, z = lon_lat_to_cartesian(
                    lon, lat, earthRadius, degrees=True
                )
                x, y, z, _, _ = subdivide_great_circle(
                    x, y, z, subdivisionResolution, earthRadius
                )
                lon, lat = cartesian_to_lon_lat(
                    x, y, z, earthRadius, degrees=True
                )
                new_coords.append([list(a) for a in zip(lon, lat)])

        if geom_type == 'LineString':
            shape = shapely.geometry.LineString(new_coords[0])
        else:
            shape = shapely.geometry.MultiLineString(new_coords)

        mask = _compute_mask_from_shapes(
            shapes1=shape,
            shapes2=polygons,
            func=_intersects,
            pool=pool,
            chunkSize=chunkSize,
            showProgress=showProgress,
        )

        masks.append(mask)
        shapes.append(shape)

    return transectNames, masks, properties, shapes


def _intersects(shape, polygons):
    tree = STRtree(polygons)
    mask = numpy.zeros(len(polygons), dtype=bool)
    indicesInShape = tree.query(shape, predicate='intersects')
    mask[indicesInShape] = True
    return mask


def _get_polygons(dsMesh, maskType):
    if maskType == 'cell':
        # polygons are vertices on cell
        vertexIndices = dsMesh.verticesOnCell.values - 1
        nVerticesOnCell = dsMesh.nEdgesOnCell.values
        maxVertices = vertexIndices.shape[1]
        for iVertex in range(maxVertices):
            mask = iVertex >= nVerticesOnCell
            # copy the last valid vertex
            vertexIndices[mask, iVertex] = vertexIndices[
                mask, nVerticesOnCell[mask] - 1
            ]
        lonVertex = dsMesh.lonVertex.values
        latVertex = dsMesh.latVertex.values
        lonCenter = dsMesh.lonCell.values
    elif maskType == 'vertex':
        # polygons are cells on vertex
        vertexIndices = dsMesh.cellsOnVertex.values - 1
        maxVertices = vertexIndices.shape[1]
        firstValid = vertexIndices[:, 0]
        for iVertex in range(1, maxVertices):
            mask = firstValid < 0
            firstValid[mask] = vertexIndices[mask, iVertex]
        assert numpy.all(firstValid >= 0)
        for iVertex in range(maxVertices):
            mask = vertexIndices[:, iVertex] < 0
            vertexIndices[mask, iVertex] = firstValid[mask]
        lonVertex = dsMesh.lonCell.values
        latVertex = dsMesh.latCell.values
        lonCenter = dsMesh.lonVertex.values
    elif maskType == 'edge':
        # oh, dear, this is a bit more complicated.  Polygons are kites
        # involving both vertices and cell centers
        verticesOnEdge = dsMesh.verticesOnEdge - 1
        cellsOnEdge = dsMesh.cellsOnEdge - 1

        nEdges = dsMesh.sizes['nEdges']
        nCells = dsMesh.sizes['nCells']
        vertexIndices = -1 * numpy.ones((nEdges, 4), int)
        vertexIndices[:, 0] = cellsOnEdge[:, 0]
        vertexIndices[:, 1] = verticesOnEdge[:, 0] + nCells
        vertexIndices[:, 2] = cellsOnEdge[:, 1]
        vertexIndices[:, 3] = verticesOnEdge[:, 1] + nCells

        # if there are invalid cells, just point to the next vertex; all
        # vertices on cell should be valid
        mask = vertexIndices[:, 0] < 0
        vertexIndices[mask, 0] = vertexIndices[mask, 1]

        mask = vertexIndices[:, 2] < 0
        vertexIndices[mask, 2] = vertexIndices[mask, 3]

        lonVertex = numpy.append(
            dsMesh.lonCell.values, dsMesh.lonVertex.values
        )
        latVertex = numpy.append(
            dsMesh.latCell.values, dsMesh.latVertex.values
        )

        lonCenter = dsMesh.lonEdge.values

    else:
        raise ValueError(f'Unknown mask type {maskType}')

    assert numpy.all(vertexIndices >= 0)

    latVertex = numpy.rad2deg(latVertex)
    # transform longitudes to [-180, 180)
    lonVertex = numpy.mod(numpy.rad2deg(lonVertex) + 180.0, 360.0) - 180.0
    lonCenter = numpy.mod(numpy.rad2deg(lonCenter) + 180.0, 360.0) - 180.0

    lon = lonVertex[vertexIndices]
    lat = latVertex[vertexIndices]

    lon, lat, duplicatePolygons = _copy_dateline_lon_lat_vertices(
        lon, lat, lonCenter
    )

    nPolygons = len(lonCenter)

    polygons = []
    for index in range(lon.shape[0]):
        coords = zip(lon[index, :], lat[index, :])
        polygons.append(shapely.geometry.Polygon(coords))

    return polygons, nPolygons, duplicatePolygons


def _copy_dateline_lon_lat_vertices(lonVertex, latVertex, lonCenter):
    nPolygons, _ = lonVertex.shape

    lonDiff = lonVertex - lonCenter.reshape(nPolygons, 1)

    # which polygons have vertices that are out of range to the west?
    outOfRange = lonDiff < -180.0
    duplicatePolygonsEast = numpy.flatnonzero(numpy.any(outOfRange, axis=1))
    lonVertex[outOfRange] += 360.0
    lonVertexToAdd = lonVertex[duplicatePolygonsEast, :] - 360.0
    latVertexToAdd = latVertex[duplicatePolygonsEast, :]

    # which polygons have vertices that are out of range to the east?
    outOfRange = lonDiff >= 180.0
    duplicatePolygonsWest = numpy.flatnonzero(numpy.any(outOfRange, axis=1))
    lonVertex[outOfRange] -= 360.0
    lonVertexToAdd = numpy.append(
        lonVertexToAdd, lonVertex[duplicatePolygonsWest, :] + 360.0, axis=0
    )
    latVertexToAdd = numpy.append(
        latVertexToAdd, latVertex[duplicatePolygonsWest, :], axis=0
    )

    lonVertex = numpy.append(lonVertex, lonVertexToAdd, axis=0)
    latVertex = numpy.append(latVertex, latVertexToAdd, axis=0)

    duplicatePolygons = numpy.append(
        duplicatePolygonsEast, duplicatePolygonsWest
    )

    # TODO: we still need to do something about cells that contain the poles

    return lonVertex, latVertex, duplicatePolygons


def _compute_seed_mask(fcSeed, lon, lat, workers):
    """
    Find the cell centers (points) closes to the given seed points and set
    the resulting mask to 1 there
    """
    points = numpy.vstack((lon, lat)).T

    tree = KDTree(points)

    mask = numpy.zeros(len(lon), dtype=numpy.int32)

    points = numpy.zeros((len(fcSeed.features), 2))
    for index, feature in enumerate(fcSeed.features):
        points[index, :] = feature['geometry']['coordinates']

    _, indices = tree.query(points, workers=workers)

    for index in indices:
        mask[index] = 1

    return mask


def _flood_fill_mask(seedMask, growMask, cellsOnCell):
    """
    Flood fill starting with a mask of seed points
    """

    maxNeighbors = cellsOnCell.shape[1]

    while True:
        neighbors = cellsOnCell[seedMask == 1, :]
        maskCount = 0
        for iNeighbor in range(maxNeighbors):
            indices = neighbors[:, iNeighbor]
            # we only want to mask valid neighbors, locations that aren't
            # already masked, and locations that we're allowed to flood
            indices = indices[indices >= 0]
            localMask = numpy.logical_and(
                seedMask[indices] == 0, growMask[indices] == 1
            )
            maskCount += numpy.count_nonzero(localMask)
            indices = indices[localMask]
            seedMask[indices] = 1

        if maskCount == 0:
            break

    return seedMask


def _compute_edge_sign(dsMesh, edgeMask, shape):
    """Compute the edge sign along a transect"""

    edge_indices = numpy.flatnonzero(edgeMask)
    voe = dsMesh.verticesOnEdge.isel(nEdges=edge_indices).values - 1

    lon = numpy.rad2deg(dsMesh.lonVertex.values[voe])
    lon = numpy.mod(lon + 180.0, 360.0) - 180.0
    lat = numpy.rad2deg(dsMesh.latVertex.values[voe])

    lonEdge = numpy.rad2deg(dsMesh.lonEdge.values[edge_indices])
    lonEdge = numpy.mod(lonEdge + 180.0, 360.0) - 180.0

    lon, lat, duplicate_edges = _copy_dateline_lon_lat_vertices(
        lon, lat, lonEdge
    )

    nEdges = dsMesh.sizes['nEdges']
    nVertices = dsMesh.sizes['nVertices']

    # give periodic copies unique edge and vertex indices
    edge_indices = numpy.append(
        edge_indices, edge_indices[duplicate_edges] + nEdges
    )

    voe = numpy.append(voe, voe[duplicate_edges, :] + nVertices, axis=0)

    unique_vertices = numpy.unique(voe.ravel())

    local_voe = numpy.zeros(voe.shape, dtype=numpy.int32)
    distance = []
    unique_lon = []
    unique_lat = []
    for local_v, v in enumerate(unique_vertices):
        local_mask = voe == v

        x = lon[local_mask][0]
        y = lat[local_mask][0]
        distance.append(shape.project(shapely.geometry.Point(x, y)))
        unique_lon.append(x)
        unique_lat.append(y)

        local_voe[local_mask] = local_v

    graph = Graph(
        n=len(unique_vertices), edges=zip(local_voe[:, 0], local_voe[:, 1])
    )
    graph.vs['distance'] = distance
    graph.vs['lon'] = unique_lon
    graph.vs['lat'] = unique_lat
    graph.vs['vertices'] = numpy.arange(len(unique_vertices))
    graph.es['edges'] = edge_indices
    graph.es['vertices'] = [(v0, v1) for v0, v1 in zip(voe[:, 0], voe[:, 1])]

    edgeSign = numpy.zeros(edgeMask.shape, dtype=numpy.int32)

    clusters = graph.connected_components()
    for cluster_index in range(len(clusters)):
        cluster = clusters.subgraph(cluster_index)
        distance = cluster.vs['distance']
        if len(cluster.es) == 1:
            edges = cluster.es.select(0)
            edge = edges[0]
            if edge.source_vertex['distance'] < edge.target_vertex['distance']:
                sign = [1]
            else:
                sign = [-1]
        else:
            start = numpy.argmin(distance)
            end = numpy.argmax(distance)
            indices = cluster.get_shortest_paths(
                v=start, to=end, output='epath'
            )[0]
            edges = cluster.es.select(indices)
            if len(edges) == 1:
                edge = edges[0]
                if (
                    edge.source_vertex['distance']
                    < edge.target_vertex['distance']
                ):
                    sign = [1]
                else:
                    sign = [-1]
            else:
                verts = numpy.array(edges['vertices'])
                sign = numpy.zeros(len(indices), dtype=numpy.int32)
                for index in range(len(indices) - 1):
                    if verts[index, 1] in verts[index + 1, :]:
                        sign[index] = 1
                    elif verts[index, 0] in verts[index + 1, :]:
                        sign[index] = -1
                    else:
                        raise ValueError(f'could not find vertex {index}')

                if verts[-1, 0] in verts[-2, :]:
                    sign[-1] = 1
                elif verts[-1, 1] in verts[-2, :]:
                    sign[-1] = -1
                else:
                    raise ValueError('could not find vertex -1')

        sign = numpy.array(sign)
        edge_indices = numpy.array(edges['edges'])
        valid = numpy.array(edge_indices) < nEdges
        edgeSign[edge_indices[valid]] = sign[valid]
        duplicate = numpy.logical_not(valid)
        edgeSign[edge_indices[duplicate] - nEdges] = sign[duplicate]

    return edgeSign
