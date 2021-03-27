import numpy
from shapely.geometry import LineString, Point


def subdivide_great_circle(x, y, z, maxRes, earthRadius):
    """
    Subdivide each segment of the transect so the horizontal resolution
    approximately matches the requested resolution

    Uses a formula for interpolating unit vectors on the sphere from
    https://en.wikipedia.org/wiki/Slerp

    Parameters
    ----------
    x, y, z : numpy.array
        The Cartesian coordinates of a transect, where the number of segments
        is ``len(x) - 1``.  ``x``, ``y`` and ``z`` are of the same length.

    maxRes : float
        The maximum allowed spacing in m after subdivision

    earthRadius : float
        The radius of the Earth in m

    Returns
    -------
    xOut, yOut, zOut : numpy.array
        The transect subdivided into segments with segment length at most
        ``maxRes``.  All the points in ``x``, ``y`` and ``z`` are guaranteed
        to be included.

    dIn : numpy.array
        The distance along the transect before subdivision

    dOut : numpy.array
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
    x, y, z : numpy.array
        Cartesian coordinates along a transect

    earth_radius : float
        The radius of the earth

    Returns
    -------
    distance : numpy.array
        The distance along the transect
    """
    distance = numpy.zeros(x.shape)
    for segIndex in range(len(x)-1):
        transectv0 = Vector(x[segIndex], y[segIndex], z[segIndex])
        transectv1 = Vector(x[segIndex+1], y[segIndex+1], z[segIndex+1])

        distance[segIndex+1] = distance[segIndex] + \
            earth_radius*angular_distance(first=transectv0, second=transectv1)

    return distance


def subdivide_planar(x, y, maxRes):
    """
    Subdivide each segment of the transect so the horizontal resolution
    approximately matches the requested resolution

    Uses a formula for interpolating unit vectors on the sphere from
    https://en.wikipedia.org/wiki/Slerp

    Parameters
    ----------
    x, y: numpy.array
        The planar coordinates of a transect, where the number of segments
        is ``len(x) - 1``.  ``x`` and ``y`` are of the same length.

    maxRes : float
        The maximum allowed spacing in m after subdivision

    Returns
    -------
    xOut, yOut : numpy.array
        The transect subdivided into segments with segment length at most
        ``maxRes``.  All the points in ``x`` and ``y`` are guaranteed to be
        included.

    dIn : numpy.array
        The distance along the transect before subdivision

    dOut : numpy.array
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


def angular_distance(x=None, y=None, z=None, first=None, second=None):
    """
    Compute angular distance between points on the sphere, following:
    https://en.wikipedia.org/wiki/Great-circle_distance

    Parameters
    ----------
    x, y, z : numpy.array, optional
        The cartsian coordinates of a transect, where the number of segments
        is ``len(x) - 1``.  ``x``, ``y`` and ``z`` are of the same length and
        all must be present if ``first`` and ``second`` are not provided.

    first, second : Vector
        The start and end points of each segment of the transect, where the
        ``x``, ``y``, and ``z`` attributes of each vector are ``numpy.array``
        objects.

    Returns
    -------
    angularDistance : numpy.array
        The angular distance (in radians) between segments of the transect.
    """
    if first is None or second is None:
        first = Vector(x[0:-1], y[0:-1], z[0:-1])
        second = Vector(x[1:], y[1:], z[1:])

    angularDistance = numpy.arctan2(_mag(_cross(first, second)),
                                    _dot(first, second))

    return angularDistance


def intersects(a1, a2, b1, b2):
    """
    Based on https://stackoverflow.com/a/26669130/7728169
    Determine if the great circle arc from ``a1`` to ``a2`` intersects that
    from ``b1`` to ``b2``.

    Parameters
    ----------
    a1, a2, b1, b2 : Vector
        Cartesian coordinates of the end points of two great circle arcs.
        The types of the attributes ``x``, ``y``, and ``z`` must either be
        ``numpy.arrays`` of identical size for all 4 vectors (in which case
        intersections are found element-wise), or scalars for
        at least one of either the ``a``s or the ``b``s.

    Returns
    -------
    intersect : numpy.array
        A boolean array of the same size as the ``a``s or the ``b``s, whichever
        is greater, indicating if the particular pair of arcs intersects
    """
    return numpy.logical_and(_straddles(a1, a2, b1, b2),
                             _straddles(b1, b2, a1, a2))


def intersection(a1, a2, b1, b2):
    """
    Based on https://stackoverflow.com/a/26669130/7728169
    Find the intersection point between great circle arc from ``a1`` to ``a2``
    and from ``b1`` to ``b2``.  The arcs should have already have been found
    to intersect by calling ``_intersects()``

    Parameters
    ----------
    a1, a2, b1, b2 : Vector
        Cartesian coordinates of the end points of two great circle arcs.
        The types of the attributes ``x``, ``y``, and ``z`` must either be
        ``numpy.arrays`` of identical size for all 4 vectors (in which case
        intersections are found element-wise), or scalars for
        at least one of either the ``a``s or the ``b``s.

    Returns
    -------
    points : Vector
        An array of Cartesian points *on the unit sphere* indicating where the
        arcs intersect
    """
    points = _cross(_cross(a1, a2), _cross(b1, b2))
    s = numpy.sign(_det(a1, b1, b2))/_mag(points)
    points = Vector(s*points.x,  s*points.y, s*points.z)
    return points


class Vector:
    """
    A class for representing Cartesian vectors with ``x``, ``y`` and ``z``
    components that are either ``float`` or ``numpy.array`` objects of identical
    size.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


def _straddles(a1, a2, b1, b2):
    """
    Based on https://stackoverflow.com/a/26669130/7728169
    Determines if the great circle segment determined by (a1, a2)
    straddles the great circle determined by (b1, b2)

    Parameters
    ----------
    a1, a2, b1, b2 : Vector
        Cartesian coordinates of the end points of two great circle arcs.
        The types of the attributes ``x``, ``y``, and ``z`` must either be
        ``numpy.arrays`` of identical size for all 4 vectors (in which case
        intersections are found element-wise), or scalars for
        at least one of either the ``a``s or the ``b``s.

    Returns
    -------
    straddle : numpy.array
        A boolean array of the same size as the ``a``s or the ``b``s, whichever
        is greater, indicating if the great circle segment determined by
        (a1, a2) straddles the great circle determined by (b1, b2)
    """
    return _det(a1, b1, b2) * _det(a2, b1, b2) < 0


def _dot(v1, v2):
    """The dot product between two ``Vector`` objects ``v1`` and ``v2``"""
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z


def _cross(v1, v2):
    """The cross product between two ``Vector`` objects ``v1`` and ``v2``"""
    return Vector(v1.y * v2.z - v1.z * v2.y,
                  v1.z * v2.x - v1.x * v2.z,
                  v1.x * v2.y - v1.y * v2.x)


def _det(v1, v2, v3):
    """The determinant of the matrix defined by the three ``Vector`` objects"""
    return _dot(v1, _cross(v2, v3))


def _mag(v):
    """The magnitude of the ``Vector`` object ``v``"""
    return numpy.sqrt(_dot(v, v))
