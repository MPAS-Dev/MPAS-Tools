from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import pyflann
from scipy import spatial
import timeit
from rasterio.features import rasterize
from affine import Affine
import shapely.geometry
from geometric_features.plot import subdivide_geom
from mpas_tools.mesh.creation.util import lonlat2xyz


def signed_distance_from_geojson(fc, lon_grd, lat_grd, earth_radius,
                                 max_length=None):
    """
    Get the distance for each point on a lon/lat grid from the closest point
    on the boundary of the geojson regions.

    Parameters
    ----------
    fc : geometrics_features.FeatureCollection
        The regions to be rasterized
    lon_grd : numpy.ndarray
        A 1D array of evenly spaced longitude values
    lat_grd : numpy.ndarray
        A 1D array of evenly spaced latitude values
    earth_radius : float
        Earth radius in meters
    max_length : float, optional
        The maximum distance (in degrees) between points on the boundary of the
        geojson region.  If the boundary is too coarse, it will be subdivided.

    Returns
    -------
    signed_distance : numpy.ndarray
       A 2D field of distances (negative inside the region, positive outside)
       to the shape boundary
    """
    distance = distance_from_geojson(fc, lon_grd, lat_grd, earth_radius,
                                     nn_search='flann', max_length=max_length)

    mask = mask_from_geojson(fc, lon_grd, lat_grd)

    signed_distance = (-2.0 * mask + 1.0) * distance
    return signed_distance


def mask_from_geojson(fc, lon_grd, lat_grd):
    """
    Make a rasterized mask on a lon/lat grid from shapes (geojson multipolygon
    data).

    Parameters
    ----------
    fc : geometrics_features.FeatureCollection
        The regions to be rasterized
    lon_grd : numpy.ndarray
        A 1D array of evenly spaced longitude values
    lat_grd : numpy.ndarray
        A 1D array of evenly spaced latitude values

    Returns
    -------
    mask : numpy.ndarray
       A 2D mask with the shapes rasterized (0.0 outside, 1.0 inside)
    """

    print("Mask from geojson")
    print("-----------------")

    nlon = len(lon_grd)
    nlat = len(lat_grd)

    dlon = (lon_grd[-1] - lon_grd[0])/(nlon-1)
    dlat = (lat_grd[-1] - lat_grd[0])/(nlat-1)

    transform = Affine(dlon, 0.0, lon_grd[0],
                       0.0, dlat, lat_grd[0])

    shapes = []
    for feature in fc.features:
        # a list of feature geometries and mask values (always 1.0)
        shape = shapely.geometry.shape(feature['geometry'])
        # expand a bit to make sure we hit the edges of the domain
        shape = shape.buffer(dlat)
        shapes.append((shapely.geometry.mapping(shape), 1.0))

    mask = rasterize(shapes, out_shape=(nlat, nlon), transform=transform)
    if np.abs(lon_grd[0] - (-180.)) < 1e-3 and \
            np.abs(lon_grd[-1] - (180.)) < 1e-3:
        # the extra column at the periodic boundary needs to be copied
        print('  fixing periodic boundary')
        mask[:, -1] = mask[:, 0]
    return mask


def distance_from_geojson(fc, lon_grd, lat_grd, earth_radius, nn_search,
                          max_length=None):
    # {{{
    """
    Get the distance for each point on a lon/lat grid from the closest point
    on the boundary of the geojson regions.

    Parameters
    ----------
    fc : geometrics_features.FeatureCollection
        The regions to be rasterized
    lon_grd : numpy.ndarray
        A 1D array of evenly spaced longitude values
    lat_grd : numpy.ndarray
        A 1D array of evenly spaced latitude values
    earth_radius : float
        Earth radius in meters
    nn_search: {'kdtree', 'flann'}
        The method used to find the nearest point on the shape boundary
    max_length : float, optional
        The maximum distance (in degrees) between points on the boundary of the
        geojson region.  If the boundary is too coarse, it will be subdivided.

    Returns
    -------
    distance : numpy.ndarray
       A 2D field of distances to the shape boundary
    """
    print("Distance from geojson")
    print("---------------------")

    print("   Finding region boundaries")
    boundary_lon = []
    boundary_lat = []
    for feature in fc.features:
        # get the boundary of each shape
        shape = shapely.geometry.shape(feature['geometry']).boundary
        if max_length is not None:
            # subdivide the shape if it's too coarse
            geom_type = shape.geom_type
            shape = subdivide_geom(shape, geom_type, max_length)
        x, y = shape.coords.xy
        boundary_lon.extend(x)
        boundary_lat.extend(y)

    boundary_lon = np.array(boundary_lon)
    boundary_lat = np.array(boundary_lat)

    # Remove point at +/- 180 lon and +/- 90 lat because these are "fake".
    # Need a little buffer (0.01 degrees) to avoid misses due to rounding.
    mask = np.logical_not(np.logical_or(
        np.logical_or(boundary_lon <= -179.99, boundary_lon >= 179.99),
        np.logical_or(boundary_lat <= -89.99, boundary_lat >= 89.99)))

    boundary_lon = boundary_lon[mask]
    boundary_lat = boundary_lat[mask]

    print("    Mean boundary latitude: {0:.2f}".format(np.mean(boundary_lat)))

    # Convert coastline points to x,y,z and create kd-tree
    npoints = len(boundary_lon)
    boundary_xyz = np.zeros((npoints, 3))
    boundary_xyz[:, 0], boundary_xyz[:, 1], boundary_xyz[:, 2] = \
        lonlat2xyz(boundary_lon, boundary_lat, earth_radius)
    flann = None
    tree = None
    if nn_search == "kdtree":
        tree = spatial.KDTree(boundary_xyz)
    elif nn_search == "flann":
        flann = pyflann.FLANN()
        flann.build_index(
            boundary_xyz,
            algorithm='kdtree',
            target_precision=1.0,
            random_seed=0)
    else:
        raise ValueError('Bad nn_search: expected kdtree or flann, got '
                         '{}'.format(nn_search))

    # Convert  backgound grid coordinates to x,y,z and put in a nx_grd x 3 array
    # for kd-tree query
    Lon_grd, Lat_grd = np.meshgrid(lon_grd, lat_grd)
    X_grd, Y_grd, Z_grd = lonlat2xyz(Lon_grd, Lat_grd, earth_radius)
    pts = np.vstack([X_grd.ravel(), Y_grd.ravel(), Z_grd.ravel()]).T

    # Find distances of background grid coordinates to the coast
    print("   Finding distance")
    start = timeit.default_timer()
    distance = None
    if nn_search == "kdtree":
        distance, _ = tree.query(pts)
    elif nn_search == "flann":
        _, distance = flann.nn_index(pts, checks=2000, random_seed=0)
        distance = np.sqrt(distance)
    end = timeit.default_timer()
    print("   Done")
    print("   {0:.0f} seconds".format(end-start))

    # Make distance array that corresponds with cell_width array
    distance = np.reshape(distance, Lon_grd.shape)

    return distance
