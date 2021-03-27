import xarray as xr
import numpy
import pyflann
import shapely.geometry
import shapely.ops
from shapely.geometry import box, Polygon, MultiPolygon, GeometryCollection
from shapely.strtree import STRtree
import progressbar
from functools import partial
import argparse

from geometric_features import read_feature_collection

from mpas_tools.transects import subdivide_great_circle, \
    lon_lat_to_cartesian, cartesian_to_lon_lat
from mpas_tools.parallel import create_pool
from mpas_tools.io import write_netcdf
from mpas_tools.logging import LoggingContext
from mpas_tools.cime.constants import constants


def compute_mpas_region_masks(dsMesh, fcMask, maskTypes=('cell', 'vertex'),
                              logger=None, pool=None, chunkSize=1000,
                              showProgress=False, subdivisionThreshold=30.):
    """
    Use shapely and processes to create a set of masks from a feature collection
    made up of regions (polygons)

    Parameters
    ----------
    dsMesh : xarray.Dataset
        An MPAS mesh on which the masks should be created

    fcMask : geometric_features.FeatureCollection
        A feature collection containing features to use to create the mask

    maskTypes : tuple of {'cell', 'edge', 'vertex'}, optional
        Which type(s) of masks to make.  Masks are created based on whether
        the latitude and longitude associated with each of these locations
        (e.g. ``dsMesh.latCell`` and ``dsMesh.lonCell`` for ``'cells'``) are
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
        lonName = 'lon{}'.format(suffix)
        latName = 'lat{}'.format(suffix)
        lat = numpy.rad2deg(dsMesh[latName].values)

        # transform longitudes to [-180, 180)
        lon = numpy.mod(numpy.rad2deg(dsMesh[lonName].values) + 180.,
                        360.) - 180.

        if logger is not None:
            logger.info('  Computing {} masks:'.format(maskType))

        # create shapely geometry for lon and lat
        points = [shapely.geometry.Point(x, y) for x, y in zip(lon, lat)]
        regionNames, masks, properties, nChar = _compute_region_masks(
            fcMask, points, logger, pool, chunkSize, showProgress,
            subdivisionThreshold)

        nPoints = len(points)

        if logger is not None:
            logger.info('  Adding masks to dataset...')
        nRegions = len(regionNames)
        # create a new data array for masks
        masksVarName = 'region{}Masks'.format(suffix)
        dsMasks[masksVarName] = \
            ((dim, 'nRegions'), numpy.zeros((nPoints, nRegions), dtype=int))

        for index in range(nRegions):
            mask = masks[index]
            dsMasks[masksVarName][:, index] = numpy.array(mask, dtype=int)

        if 'regionNames' not in dsMasks:
            # create a new data array for mask names
            dsMasks['regionNames'] = (('nRegions',),
                                      numpy.zeros((nRegions,),
                                                  dtype='|S{}'.format(nChar)))

            for index in range(nRegions):
                dsMasks['regionNames'][index] = regionNames[index]

        for propertyName in properties:
            if propertyName not in dsMasks:
                dsMasks[propertyName] = (('nRegions',),
                                         properties[propertyName])
    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def entry_point_compute_mpas_region_masks():
    """ Entry point for ``compute_mpas_region_masks()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mesh_file_name", dest="mesh_file_name",
                        type=str, required=True,
                        help="An MPAS mesh file")
    parser.add_argument("-g", "--geojson_file_name",
                        dest="geojson_file_name", type=str, required=True,
                        help="An Geojson file containing mask regions")
    parser.add_argument("-o", "--mask_file_name", dest="mask_file_name",
                        type=str, required=True,
                        help="An output MPAS region masks file")
    parser.add_argument("-t", "--mask_types", nargs='+', dest="mask_types",
                        type=str,
                        help="Which type(s) of masks to make: cell, edge or "
                             "vertex.  Default is cell and vertex.")
    parser.add_argument("-c", "--chunk_size", dest="chunk_size", type=int,
                        default=1000,
                        help="The number of cells, vertices or edges that are "
                             "processed in one operation")
    parser.add_argument("--show_progress", dest="show_progress",
                        action="store_true",
                        help="Whether to show a progress bar")
    parser.add_argument("-s", "--subdivision", dest="subdivision", type=float,
                        default=30.,
                        help="A threshold in degrees (lon or lat) above which "
                             "the mask region will be subdivided into smaller "
                             "polygons for faster intersection checking")
    parser.add_argument(
        "--process_count", required=False, dest="process_count", type=int,
        help="The number of processes to use to compute masks.  The "
             "default is to use all available cores")
    parser.add_argument(
        "--multiprocessing_method", dest="multiprocessing_method",
        default='forkserver',
        help="The multiprocessing method use for python mask creation "
             "('fork', 'spawn' or 'forkserver')")
    args = parser.parse_args()

    dsMesh = xr.open_dataset(args.mesh_file_name, decode_cf=False,
                             decode_times=False)
    fcMask = read_feature_collection(args.geojson_file_name)

    pool = create_pool(process_count=args.process_count,
                       method=args.multiprocessing_method)

    if args.mask_types is None:
        args.mask_types = ('cell', 'vertex')

    with LoggingContext('compute_mpas_region_masks') as logger:
        dsMasks = compute_mpas_region_masks(
            dsMesh=dsMesh, fcMask=fcMask, maskTypes=args.mask_types,
            logger=logger, pool=pool, chunkSize=args.chunk_size,
            showProgress=args.show_progress,
            subdivisionThreshold=args.subdivision)

    write_netcdf(dsMasks, args.mask_file_name)


def compute_mpas_transect_masks(dsMesh, fcMask, earthRadius,
                                maskTypes=('cell', 'edge', 'vertex'),
                                logger=None, pool=None, chunkSize=1000,
                                showProgress=False, subdivisionResolution=10e3):
    """
    Use shapely and processes to create a set of masks from a feature collection
    made up of transects (line strings)

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
        (e.g. ``dsMesh.latCell`` and ``dsMesh.lonCell`` for ``'cells'``) are
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
            logger.info('  Computing {} masks:'.format(maskType))

        polygons, nPolygons, duplicatePolygons = _get_polygons(dsMesh, maskType)
        transectNames, masks, properties, nChar = _compute_transect_masks(
            fcMask, polygons, logger, pool, chunkSize, showProgress,
            subdivisionResolution, earthRadius)

        if logger is not None:
            logger.info('  Adding masks to dataset...')
        nTransects = len(transectNames)
        # create a new data array for masks
        masksVarName = 'transect{}Masks'.format(suffix)
        dsMasks[masksVarName] = \
            ((dim, 'nTransects'),
             numpy.zeros((nPolygons, nTransects), dtype=int))

        for index in range(nTransects):
            maskAndDuplicates = masks[index]
            mask = maskAndDuplicates[0:nPolygons]

            mask[duplicatePolygons] = \
                numpy.logical_or(mask[duplicatePolygons],
                                 maskAndDuplicates[nPolygons:])
            dsMasks[masksVarName][:, index] = numpy.array(mask, dtype=int)

        if 'transectNames' not in dsMasks:
            # create a new data array for mask names
            dsMasks['transectNames'] = (('nTransects',),
                                        numpy.zeros((nTransects,),
                                                    dtype='|S{}'.format(nChar)))

            for index in range(nTransects):
                dsMasks['transectNames'][index] = transectNames[index]

        for propertyName in properties:
            if propertyName not in dsMasks:
                dsMasks[propertyName] = (('nTransects',),
                                         properties[propertyName])
    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def entry_point_compute_mpas_transect_masks():
    """ Entry point for ``compute_mpas_transect_masks()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mesh_file_name", dest="mesh_file_name",
                        type=str, required=True,
                        help="An MPAS mesh file")
    parser.add_argument("-g", "--geojson_file_name",
                        dest="geojson_file_name", type=str, required=True,
                        help="An Geojson file containing transects")
    parser.add_argument("-o", "--mask_file_name", dest="mask_file_name",
                        type=str, required=True,
                        help="An output MPAS transect masks file")
    parser.add_argument("-t", "--mask_types", nargs='+', dest="mask_types",
                        type=str,
                        help="Which type(s) of masks to make: cell, edge or "
                             "vertex.  Default is cell, edge and vertex.")
    parser.add_argument("-c", "--chunk_size", dest="chunk_size", type=int,
                        default=1000,
                        help="The number of cells, vertices or edges that are "
                             "processed in one operation")
    parser.add_argument("--show_progress", dest="show_progress",
                        action="store_true",
                        help="Whether to show a progress bar")
    parser.add_argument("-s", "--subdivision", dest="subdivision", type=float,
                        help="The maximum resolution (in meters) of segments "
                             "in a transect.  If a transect is too coarse, it "
                             "will be subdivided.  Default is no subdivision.")
    parser.add_argument(
        "--process_count", required=False, dest="process_count", type=int,
        help="The number of processes to use to compute masks.  The "
             "default is to use all available cores")
    parser.add_argument(
        "--multiprocessing_method", dest="multiprocessing_method",
        default='forkserver',
        help="The multiprocessing method use for python mask creation "
             "('fork', 'spawn' or 'forkserver')")
    args = parser.parse_args()

    dsMesh = xr.open_dataset(args.mesh_file_name, decode_cf=False,
                             decode_times=False)
    fcMask = read_feature_collection(args.geojson_file_name)

    pool = create_pool(process_count=args.process_count,
                       method=args.multiprocessing_method)

    if args.mask_types is None:
        args.mask_types = ('cell', 'edge', 'vertex')

    earth_radius = constants['SHR_CONST_REARTH']

    with LoggingContext('compute_mpas_transect_masks') as logger:
        dsMasks = compute_mpas_transect_masks(
            dsMesh=dsMesh, fcMask=fcMask, earthRadius=earth_radius,
            maskTypes=args.mask_types, logger=logger, pool=pool,
            chunkSize=args.chunk_size, showProgress=args.show_progress,
            subdivisionResolution=args.subdivision)

    write_netcdf(dsMasks, args.mask_file_name)


def compute_mpas_flood_fill_mask(dsMesh, fcSeed, logger=None):
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

    logger : logging.Logger, optional
        A logger for the output if not stdout

    Returns
    -------
    dsMask : xarray.Dataset
        The masks

    """

    dsMasks = xr.Dataset()

    lat = numpy.rad2deg(dsMesh.latCell.values)

    # transform longitudes to [-180, 180)
    lon = numpy.mod(numpy.rad2deg(dsMesh.lonCell.values) + 180.,
                    360.) - 180.

    if logger is not None:
        logger.info('  Computing flood fill mask on cells:')

    mask = _compute_seed_mask(fcSeed, lon, lat)

    cellsOnCell = dsMesh.cellsOnCell.values - 1

    mask = _flood_fill_mask(mask, cellsOnCell)

    if logger is not None:
        logger.info('  Adding masks to dataset...')
    # create a new data array for the mask
    masksVarName = 'cellSeedMask'
    dsMasks[masksVarName] = (('nCells',), numpy.array(mask, dtype=int))

    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def entry_point_compute_mpas_flood_fill_mask():
    """ Entry point for ``compute_mpas_flood_fill_mask()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mesh_file_name", dest="mesh_file_name",
                        type=str, required=True,
                        help="An MPAS mesh file")
    parser.add_argument("-g", "--geojson_file_name",
                        dest="geojson_file_name", type=str, required=True,
                        help="An Geojson file containing points at which to "
                             "start the flood fill")
    parser.add_argument("-o", "--mask_file_name", dest="mask_file_name",
                        type=str, required=True,
                        help="An output MPAS region masks file")
    args = parser.parse_args()

    dsMesh = xr.open_dataset(args.mesh_file_name, decode_cf=False,
                             decode_times=False)
    fcSeed = read_feature_collection(args.geojson_file_name)

    with LoggingContext('compute_mpas_flood_fill_mask') as logger:
        dsMasks = compute_mpas_flood_fill_mask(
            dsMesh=dsMesh, fcSeed=fcSeed, logger=logger)

    write_netcdf(dsMasks, args.mask_file_name)


def compute_lon_lat_region_masks(lon, lat, fcMask, logger=None, pool=None,
                                 chunkSize=1000, showProgress=False,
                                 subdivisionThreshold=30.):
    """
    Use shapely and processes to create a set of masks from a feature collection
    made up of regions (polygons) on a tensor lon/lat grid

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

    Lon, Lat = numpy.meshgrid(lon, lat)

    shape = Lon.shape

    Lon = Lon.ravel()
    Lat = Lat.ravel()

    # create shapely geometry for lon and lat
    points = [shapely.geometry.Point(x, y) for x, y in zip(Lon, Lat)]
    regionNames, masks, properties, nChar = _compute_region_masks(
        fcMask, points, logger, pool, chunkSize, showProgress,
        subdivisionThreshold)

    nlon = len(lon)
    nlat = len(lat)

    if logger is not None:
        logger.info('  Adding masks to dataset...')
    nRegions = len(regionNames)
    # create a new data array for masks
    masksVarName = 'regionMasks'
    dsMasks[masksVarName] = \
        (('lat', 'lon', 'nRegions'), numpy.zeros((nlat, nlon, nRegions),
                                                 dtype=int))

    for index in range(nRegions):
        mask = masks[index]
        dsMasks[masksVarName][:, :, index] = \
            numpy.array(mask.reshape(shape), dtype=int)

    # create a new data array for mask names
    dsMasks['regionNames'] = (('nRegions',),
                              numpy.zeros((nRegions,),
                                          dtype='|S{}'.format(nChar)))

    for index in range(nRegions):
        dsMasks['regionNames'][index] = regionNames[index]

    for propertyName in properties:
        if propertyName not in dsMasks:
            dsMasks[propertyName] = (('nRegions',),
                                     properties[propertyName])
    if logger is not None:
        logger.info('  Done.')

    return dsMasks


def _get_region_names_and_properties(fc):
    regionNames = []
    for feature in fc.features:
        name = feature['properties']['name']
        regionNames.append(name)

    propertyNames = set()
    for feature in fc.features:
        for propertyName in feature['properties']:
            if propertyName not in ['name', 'author', 'tags', 'component',
                                    'object']:
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


def entry_point_compute_lon_lat_region_masks():
    """ Entry point for ``compute_lon_lat_region_masks()``"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--grid_file_name", dest="grid_file_name",
                        type=str, required=True,
                        help="An input lon/lat grid file")
    parser.add_argument("--lon", dest="lon", default="lon", type=str,
                        help="The name of the longitude coordinate")
    parser.add_argument("--lat", dest="lat", default="lat", type=str,
                        help="The name of the latitude coordinate")
    parser.add_argument("-g", "--geojson_file_name",
                        dest="geojson_file_name", type=str, required=True,
                        help="An Geojson file containing mask regions")
    parser.add_argument("-o", "--mask_file_name", dest="mask_file_name",
                        type=str, required=True,
                        help="An output MPAS region masks file")
    parser.add_argument("-c", "--chunk_size", dest="chunk_size", type=int,
                        default=1000,
                        help="The number of grid points that are "
                             "processed in one operation")
    parser.add_argument("--show_progress", dest="show_progress",
                        action="store_true",
                        help="Whether to show a progress bar")
    parser.add_argument("-s", "--subdivision", dest="subdivision", type=float,
                        default=30.,
                        help="A threshold in degrees (lon or lat) above which "
                             "the mask region will be subdivided into smaller "
                             "polygons for faster intersection checking")
    parser.add_argument(
        "--process_count", required=False, dest="process_count", type=int,
        help="The number of processes to use to compute masks.  The "
             "default is to use all available cores")
    parser.add_argument(
        "--multiprocessing_method", dest="multiprocessing_method",
        default='forkserver',
        help="The multiprocessing method use for python mask creation "
             "('fork', 'spawn' or 'forkserver')")
    args = parser.parse_args()

    dsGrid = xr.open_dataset(args.grid_file_name, decode_cf=False,
                             decode_times=False)
    lon = dsGrid[args.lon].values
    lat = dsGrid[args.lat].values

    fcMask = read_feature_collection(args.geojson_file_name)

    pool = create_pool(process_count=args.process_count,
                       method=args.multiprocessing_method)

    with LoggingContext('compute_lon_lat_region_masks') as logger:
        dsMasks = compute_lon_lat_region_masks(
            lon=lon, lat=lat, fcMask=fcMask, logger=logger, pool=pool,
            chunkSize=args.chunk_size, showProgress=args.show_progress,
            subdivisionThreshold=args.subdivision)

    write_netcdf(dsMasks, args.mask_file_name)


def _compute_mask_from_shapes(shapes1, shapes2, func, pool, chunkSize,
                              showProgress):
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
            widgets = ['    ', progressbar.Percentage(), ' ',
                       progressbar.Bar(), ' ', progressbar.ETA()]
            bar = progressbar.ProgressBar(widgets=widgets,
                                          maxval=nChunks).start()
        else:
            bar = None

        mask = numpy.zeros((nShapes2,), bool)
        for iChunk, maskChunk in \
                enumerate(pool.imap(partial_func, chunks)):
            mask[indices[iChunk]:indices[iChunk + 1]] = maskChunk
            if showProgress:
                bar.update(iChunk + 1)
        if showProgress:
            bar.finish()
    return mask


def _compute_region_masks(fcMask, points, logger, pool, chunkSize, showProgress,
                          threshold):
    """
    Build a region mask file from the given mesh and geojson file defining
    a set of regions.
    """

    regionNames, properties = _get_region_names_and_properties(fcMask)

    masks = []

    nChar = 0
    for feature in fcMask.features:
        name = feature['properties']['name']

        if logger is not None:
            logger.info('    {}'.format(name))

        shape = shapely.geometry.shape(feature['geometry'])
        shapes = _katana(shape, threshold=threshold)

        mask = _compute_mask_from_shapes(
            shapes1=shapes, shapes2=points, func=_contains,
            pool=pool, chunkSize=chunkSize, showProgress=showProgress)

        nChar = max(nChar, len(name))

        masks.append(mask)

    return regionNames, masks, properties, nChar


def _contains(shapes, points):
    tree = STRtree(points)
    indices = dict((id(point), index) for index, point in enumerate(points))
    mask = numpy.zeros(len(points), dtype=bool)
    for shape in shapes:
        pointsToCheck = tree.query(shape)
        indicesInShape = [indices[id(point)] for point in pointsToCheck if
                          shape.contains(point)]
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
    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

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
        a = box(bounds[0], bounds[1], bounds[2], bounds[1]+height/2)
        b = box(bounds[0], bounds[1]+height/2, bounds[2], bounds[3])
    else:
        # split top to bottom
        a = box(bounds[0], bounds[1], bounds[0]+width/2, bounds[3])
        b = box(bounds[0]+width/2, bounds[1], bounds[2], bounds[3])
    result = []
    for d in (a, b,):
        c = geometry.intersection(d)
        if not isinstance(c, GeometryCollection):
            c = [c]
        for e in c:
            if isinstance(e, (Polygon, MultiPolygon)):
                result.extend(_katana(e, threshold, count+1, maxcount))
    if count > 0:
        return result
    # convert multipart into singlepart
    final_result = []
    for g in result:
        if isinstance(g, MultiPolygon):
            final_result.extend(g)
        else:
            final_result.append(g)
    return final_result


def _compute_transect_masks(fcMask, polygons, logger, pool, chunkSize,
                            showProgress, subdivisionResolution, earthRadius):
    """
    Build a transect mask file from the given mesh and geojson file defining
    a set of transects.
    """

    transectNames, properties = _get_region_names_and_properties(fcMask)

    masks = []

    nChar = 0
    for feature in fcMask.features:
        name = feature['properties']['name']

        if logger is not None:
            logger.info('    {}'.format(name))

        geom_type = feature['geometry']['type']
        if geom_type == 'LineString':
            coordinates = [feature['geometry']['coordinates']]
        elif geom_type == 'MultiLineString':
            coordinates = feature['geometry']['coordinates']
        else:
            raise ValueError('Unexpected geometry type {}'.format(geom_type))

        new_coords = []
        for coord_index, coords in enumerate(coordinates):

            if subdivisionResolution is None:
                new_coords.append(coords)
            else:
                lon, lat = zip(*coords)
                x, y, z = lon_lat_to_cartesian(
                    lon, lat, earthRadius, degrees=True)
                x, y, z, _, _ = subdivide_great_circle(
                    x, y, z, subdivisionResolution, earthRadius)
                lon, lat = cartesian_to_lon_lat(
                    x, y, z, earthRadius, degrees=True)
                new_coords.append([list(a) for a in zip(lon, lat)])

        if geom_type == 'LineString':
            shape = shapely.geometry.LineString(new_coords[0])
        else:
            shape = shapely.geometry.MultiLineString(new_coords)

        mask = _compute_mask_from_shapes(
            shapes1=shape, shapes2=polygons, func=_intersects,
            pool=pool, chunkSize=chunkSize, showProgress=showProgress)

        nChar = max(nChar, len(name))

        masks.append(mask)

    return transectNames, masks, properties, nChar


def _intersects(shape, polygons):
    tree = STRtree(polygons)
    indices = dict((id(polygon), index) for index, polygon in
                   enumerate(polygons))
    mask = numpy.zeros(len(polygons), dtype=bool)
    polygonsToCheck = tree.query(shape)
    indicesInShape = [indices[id(polygon)] for polygon in polygonsToCheck if
                      shape.intersects(polygon)]
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
            vertexIndices[mask, iVertex] = \
                vertexIndices[mask, nVerticesOnCell[mask]-1]
        lonVertex = dsMesh.lonVertex.values
        latVertex = dsMesh.latVertex.values
        lonCenter = dsMesh.lonCell.values
    elif maskType == 'vertex':
        # polygons are cells on vertex
        vertexIndices = dsMesh.cellsOnVertex.values - 1
        valid = numpy.amax(vertexIndices >= 0, axis=1)
        vertexIndices = vertexIndices[valid, :]
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

        lonVertex = numpy.append(dsMesh.lonCell.values,
                                 dsMesh.lonVertex.values)
        latVertex = numpy.append(dsMesh.latCell.values,
                                 dsMesh.latVertex.values)

        lonCenter = dsMesh.lonEdge.values

    else:
        raise ValueError('Unknown mask type {}'.format(maskType))

    assert numpy.all(vertexIndices >= 0)

    latVertex = numpy.rad2deg(latVertex)
    # transform longitudes to [-180, 180)
    lonVertex = numpy.mod(numpy.rad2deg(lonVertex) + 180., 360.) - 180.
    lonCenter = numpy.mod(numpy.rad2deg(lonCenter) + 180., 360.) - 180.

    lon = lonVertex[vertexIndices]
    lat = latVertex[vertexIndices]

    lon, lat, duplicatePolygons = _copy_dateline_lon_lat_vertices(lon, lat,
                                                                  lonCenter)

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
    outOfRange = lonDiff < -180.
    duplicatePolygonsEast = numpy.any(outOfRange, axis=1)
    lonVertex[outOfRange] += 360.
    lonVertexToAdd = lonVertex[duplicatePolygonsEast, :] - 360.
    latVertexToAdd = latVertex[duplicatePolygonsEast, :]

    # which polygons have vertices that are out of range to the east?
    outOfRange = lonDiff >= 180.
    duplicatePolygonsWest = numpy.any(outOfRange, axis=1)
    lonVertex[outOfRange] -= 360.
    lonVertexToAdd = numpy.append(lonVertexToAdd,
                                  lonVertex[duplicatePolygonsWest, :] + 360.,
                                  axis=0)
    latVertexToAdd = numpy.append(latVertexToAdd,
                                  latVertex[duplicatePolygonsWest, :],
                                  axis=0)

    lonVertex = numpy.append(lonVertex, lonVertexToAdd, axis=0)
    latVertex = numpy.append(latVertex, latVertexToAdd, axis=0)

    duplicatePolygons = numpy.logical_or(duplicatePolygonsEast,
                                         duplicatePolygonsWest)

    # TODO: we still need to do something about cells that contain the poles

    return lonVertex, latVertex, duplicatePolygons


def _compute_seed_mask(fcSeed, lon, lat):
    """
    Find the cell centers (points) closes to the given seed points and set
    the resulting mask to 1 there
    """
    points = numpy.vstack((lon, lat)).T
    flann = pyflann.FLANN()
    flann.build_index(points, algorithm='kmeans', target_precision=1.0,
                      random_seed=0)

    mask = numpy.zeros(len(lon), dtype=int)

    points = numpy.zeros((len(fcSeed.features), 2))
    for index, feature in enumerate(fcSeed.features):
        points[index, :] = feature['geometry']['coordinates']

    indices, distances = flann.nn_index(points, checks=2000, random_seed=0)
    for index in indices:
        mask[index] = 1

    return mask


def _flood_fill_mask(mask, cellsOnCell):
    """
    Flood fill starting with a mask of seed points
    """

    maxNeighbors = cellsOnCell.shape[1]

    while True:
        neighbors = cellsOnCell[mask == 1, :]
        maskCount = 0
        for iNeighbor in range(maxNeighbors):
            indices = neighbors[:, iNeighbor]
            # we only want to mask valid neighbors and locations that aren't
            # already masked
            indices = indices[indices >= 0]
            localMask = mask[indices] == 0
            maskCount += numpy.count_nonzero(localMask)
            indices = indices[localMask]
            mask[indices] = 1

        if maskCount == 0:
            break

    return mask
