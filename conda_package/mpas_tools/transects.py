import numpy
from shapely.geometry import LineString, Point

from mpas_tools.vector import Vector


def subdivide_great_circle(x, y, z, maxRes, earthRadius):
    """
    Subdivide each segment of the transect so the horizontal resolution
    approximately matches the requested resolution

    Uses a formula for interpolating unit vectors on the sphere from
    https://en.wikipedia.org/wiki/Slerp

    Parameters
    ----------
    x : numpy.ndarray
        The Cartesian x coordinate of a transect, where the number of segments
        is ``len(x) - 1``.  ``x``, ``y`` and ``z`` are of the same length.
    y : numpy.ndarray
        The Cartesian y coordinate of the transect
    z : numpy.ndarray
        The Cartesian z coordinate of the transect

    maxRes : float
        The maximum allowed spacing in m after subdivision

    earthRadius : float
        The radius of the Earth in m

    Returns
    -------
    xOut : numpy.ndarray
        The Cartesian x values of the transect subdivided into segments with
        segment length at most ``maxRes``.  All the points in ``x``, ``y`` and
        ``z`` are guaranteed to be included.

    yOut : numpy.ndarray
        The Cartesian y values of the subdivided transect

    zOut : numpy.ndarray
        The Cartesian y values of the subdivided transect

    dIn : numpy.ndarray
        The distance along the transect before subdivision

    dOut : numpy.ndarray
        The distance along the transect after subdivision

    """

    angularDistance = angular_distance(x=x, y=y, z=z)

    dx = angularDistance * earthRadius

    nSegments = numpy.maximum(
        (dx / maxRes + 0.5).astype(int), 1)

    dIn = numpy.zeros(x.shape)
    dIn[1:] = numpy.cumsum(dx)

    frac = []
    outIndices = []
    delta = []
    for index in range(len(dIn) - 1):
        n = nSegments[index]
        frac.extend(numpy.arange(0, n)/n)
        outIndices.extend(index*numpy.ones(n, int))
        delta.extend(angularDistance[index]*numpy.ones(n))
    frac.append(1.)
    outIndices.append(len(dIn) - 2)
    delta.append(angularDistance[-1])

    frac = numpy.array(frac)
    delta = numpy.array(delta)
    outIndices = numpy.array(outIndices)

    a = numpy.ones(delta.shape)
    b = numpy.zeros(delta.shape)
    mask = delta > 0.
    denom = 1./numpy.sin(delta[mask])
    a[mask] = denom*numpy.sin((1.-frac[mask])*delta[mask])
    b[mask] = denom*numpy.sin(frac[mask]*delta[mask])

    xOut = a*x[outIndices] + b*x[outIndices+1]
    yOut = a*y[outIndices] + b*y[outIndices+1]
    zOut = a*z[outIndices] + b*z[outIndices+1]

    dOut = (-frac + 1.)*dIn[outIndices] + frac*dIn[outIndices+1]

    return xOut, yOut, zOut, dIn, dOut


def cartesian_to_great_circle_distance(x, y, z, earth_radius):
    """
    Cartesian transect points to great-circle distance

    Parameters
    ----------
    x : numpy.ndarray
        The Cartesian x coordinate of a transect
    y : numpy.ndarray
        The Cartesian y coordinate of the transect
    z : numpy.ndarray
        The Cartesian z coordinate of the transect

    earth_radius : float
        The radius of the earth

    Returns
    -------
    distance : numpy.ndarray
        The distance along the transect
    """
    distance = numpy.zeros(x.shape)
    for segIndex in range(len(x)-1):
        transectv0 = Vector(x[segIndex], y[segIndex], z[segIndex])
        transectv1 = Vector(x[segIndex+1], y[segIndex+1], z[segIndex+1])

        distance[segIndex+1] = distance[segIndex] + \
            earth_radius*transectv0.angular_distance(transectv1)

    return distance


def subdivide_planar(x, y, maxRes):
    """
    Subdivide each segment of the transect so the horizontal resolution
    approximately matches the requested resolution

    Uses a formula for interpolating unit vectors on the sphere from
    https://en.wikipedia.org/wiki/Slerp

    Parameters
    ----------
    x : numpy.ndarray
        The planar x coordinate of a transect, where the number of segments
        is ``len(x) - 1``

    y : numpy.ndarray
        The planar y coordinates of the transect, the same length as ``x``

    maxRes : float
        The maximum allowed spacing in m after subdivision

    Returns
    -------
    xOut : numpy.ndarray
        The x coordinate of the transect, subdivided into segments with length
        at most ``maxRes``.  All the points in ``x`` are guaranteed to be
        included.

    yOut : numpy.ndarray
        The y coordinate of the transect.  All the points in ``y`` are
        guaranteed to be included.

    dIn : numpy.ndarray
        The distance along the transect before subdivision

    dOut : numpy.ndarray
        The distance along the transect after subdivision
    """

    dx = numpy.zeros(len(x)-1)
    for index in range(len(x)-1):
        start = Point(x[index], y[index])
        end = Point(x[index+1], y[index+1])
        segment = LineString([start, end])
        dx[index] = segment.length

    nSegments = numpy.maximum(
        (dx / maxRes + 0.5).astype(int), 1)

    dIn = numpy.zeros(x.shape)
    dIn[1:] = numpy.cumsum(dx)

    frac = []
    outIndices = []
    for index in range(len(dIn) - 1):
        n = nSegments[index]
        frac.extend(numpy.arange(0, n)/n)
        outIndices.extend(index*numpy.ones(n, int))
    frac.append(1.)
    outIndices.append(len(dIn) - 2)

    frac = numpy.array(frac)
    outIndices = numpy.array(outIndices)

    xOut = (-frac + 1.)*x[outIndices] + frac*x[outIndices+1]
    yOut = (-frac + 1.)*y[outIndices] + frac*y[outIndices+1]
    dOut = (-frac + 1.)*dIn[outIndices] + frac*dIn[outIndices+1]

    return xOut, yOut, dIn, dOut


def lon_lat_to_cartesian(lon, lat, earth_radius, degrees):
    """Convert from lon/lat to Cartesian x, y, z"""

    if degrees:
        lon = numpy.deg2rad(lon)
        lat = numpy.deg2rad(lat)
    x = earth_radius * numpy.cos(lat) * numpy.cos(lon)
    y = earth_radius * numpy.cos(lat) * numpy.sin(lon)
    z = earth_radius * numpy.sin(lat)
    return x, y, z


def cartesian_to_lon_lat(x, y, z, earth_radius, degrees):
    """Convert from  Cartesian x, y, z to lon/lat"""
    lon = numpy.arctan2(y, x)
    lat = numpy.arcsin(z/earth_radius)
    if degrees:
        lon = numpy.rad2deg(lon)
        lat = numpy.rad2deg(lat)
    return lon, lat


def angular_distance(x, y, z):
    """
    Compute angular distance between points on the sphere, following:
    https://en.wikipedia.org/wiki/Great-circle_distance

    Parameters
    ----------
    x : numpy.ndarray
        The Cartesian x coordinate of a transect, where the number of segments
        is ``len(x) - 1``.  ``x``, ``y`` and ``z`` are of the same length.

    y : numpy.ndarray
        The Cartesian y coordinate of the transect

    z : numpy.ndarray
        The Cartesian z coordinate of the transect

    Returns
    -------
    distance : numpy.ndarray
        The angular distance (in radians) between segments of the transect.
    """
    first = Vector(x[0:-1], y[0:-1], z[0:-1])
    second = Vector(x[1:], y[1:], z[1:])

    distance = first.angular_distance(second)
    return distance
