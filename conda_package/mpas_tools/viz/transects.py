import numpy
import xarray

from scipy.spatial import cKDTree


def make_triangle_tree(dsTris):
    """
    Make a KD-Tree for finding triangle edges that are near enough to transect
    segments that they might intersect

    Parameters
    ----------
    dsTris : xarray.Dataset
        A dataset that defines triangles, the results of calling
        ``mesh_to_triangles()``

    Returns
    -------
    tree : scipy.spatial.cKDTree
        A tree of edge centers from triangles making up an MPAS mesh
    """

    nTriangles = dsTris.sizes['nTriangles']
    nNodes = dsTris.sizes['nNodes']
    nodeCoords = numpy.zeros((nTriangles*nNodes, 3))
    nodeCoords[:, 0] = dsTris.xNode.values.ravel()
    nodeCoords[:, 1] = dsTris.yNode.values.ravel()
    nodeCoords[:, 2] = dsTris.zNode.values.ravel()

    nextTri, nextNode = numpy.meshgrid(
        numpy.arange(nTriangles), numpy.mod(numpy.arange(nNodes) + 1, 3),
        indexing='ij')
    nextIndices = nNodes*nextTri.ravel() + nextNode.ravel()

    # edge centers are half way between adjacent nodes (ignoring great-circle
    # distance)
    edgeCoords = 0.5*(nodeCoords + nodeCoords[nextIndices, :])

    tree = cKDTree(data=edgeCoords, copy_data=True)
    return tree


def find_transect_cells_and_weights(lonTransect, latTransect, dsTris, dsMesh,
                                    tree, degrees=True, subdivisionRes=10e3):
    """
    Find "nodes" where the transect intersects the edges of the triangles
    that make up MPAS cells.

    Parameters
    ----------
    lonTransect, latTransect : xarray.DataArray
        The latitude and longitude of segments making up the transect

    dsTris : xarray.Dataset
        A dataset that defines triangles, the results of calling
        ``mesh_to_triangles()``

    dsMesh : xarray.Dataset
        A data set with the full MPAS mesh.

    tree : scipy.spatial.cKDTree
        A tree of edge centers from triangles making up an MPAS mesh, the return
        value from ``make_triangle_tree()``

    degrees : bool, optional
        Whether ``lonTransect`` and ``latTransect`` are in degrees (as opposed
        to radians).

    subdivisionRes : float, optional
        Resolution in m to use to subdivide the transect when looking for
        intersection candidates.  Should be small enough that curvature is
        small.

    Returns
    -------
    dsOut : xarray.Dataset
        A dataset that contains "nodes" where the transect intersects the
        edges of the triangles in ``dsTris``.  The nodes also includes the two
        end points of the transect, which typically lie within triangles. Each
        internal node (that is, not including the end points) is purposefully
        repeated twice, once for each triangle that node touches.  This allows
        for discontinuous fields between triangles (e.g. if one wishes to plot
        constant values on each MPAS cell).  The Cartesian and lon/lat
        coordinates of these nodes are ``xCartNode``, ``yCartNode``, `
        `zCartNode``, ``lonNode`` and ``latNode``.  The distance along the
        transect of each intersection is ``dNode``. The index of the triangle
        and the first triangle node in ``dsTris`` associated with each
        intersection node are given by ``horizTriangleIndices`` and
        ``horizTriangleNodeIndices``, respectively. The second node on the
        triangle for the edge associated with the intersection is given by
        ``numpy.mod(horizTriangleNodeIndices + 1, 3)``.

        The MPAS cell tha a given node belongs to is given by ``cellIndices``.
        Each node also has an associated set of 6 ``interpHorizCellIndices`` and
        ``interpHorizCellWeights`` that can be used to interpolate from MPAS
        cell centers to nodes first with area-weighted averaging to MPAS
        vertices and then linear interpolation along triangle edges.  Some of
        the weights may be zero, in which case the associated
        ``interpHorizCellIndices`` will be -1.

        Finally, ``lonTransect`` and ``latTransect`` are included in the
        dataset, along with Cartesian coordinates ``xCartTransect``,
        ``yCartTransect``, `zCartTransect``, and ``dTransect``, the great-circle
        distance along the transect of each original transect point.  In order
        to interpolate values (e.g. observations) from the original transect
        points to the intersection nodes, linear interpolation indices
        ``transectIndicesOnHorizNode`` and weights
        ``transectWeightsOnHorizNode`` are provided.  The values at nodes are
        found by::

          nodeValues = ((transectValues[transectIndicesOnHorizNode] *
                         transectWeightsOnHorizNode)
                        + (transectValues[transectIndicesOnHorizNode+1] *
                           (1.0 - transectWeightsOnHorizNode))
    """
    earth_radius = dsMesh.attrs['sphere_radius']
    buffer = numpy.maximum(numpy.amax(dsMesh.dvEdge.values),
                           numpy.amax(dsMesh.dcEdge.values))

    x, y, z = _lon_lat_to_cartesian(lonTransect, latTransect, earth_radius,
                                    degrees)

    nTriangles = dsTris.sizes['nTriangles']
    nNodes = dsTris.sizes['nNodes']
    nodeCellWeights = dsTris.nodeCellWeights.values
    nodeCellIndices = dsTris.nodeCellIndices.values

    xNode = dsTris.xNode.values.ravel()
    yNode = dsTris.yNode.values.ravel()
    zNode = dsTris.zNode.values.ravel()

    dTransect = numpy.zeros(lonTransect.shape)

    dNode = None
    xOut = None
    yOut = None
    zOut = None
    tris = None
    nodes = None
    interpCells = None
    cellWeights = None

    nHorizWeights = 6

    first = True

    dStart = 0.
    for segIndex in range(len(x)-1):
        transectv0 = Vector(x[segIndex].values,
                            y[segIndex].values,
                            z[segIndex].values)
        transectv1 = Vector(x[segIndex+1].values,
                            y[segIndex+1].values,
                            z[segIndex+1].values)

        subSlice = slice(segIndex, segIndex+2)
        xSub, ySub, zSub, _, _ = subdivide_great_circle(
            x[subSlice].values, y[subSlice].values, z[subSlice].values,
            subdivisionRes, earth_radius)

        coords = numpy.zeros((len(xSub), 3))
        coords[:, 0] = xSub
        coords[:, 1] = ySub
        coords[:, 2] = zSub
        radius = buffer + subdivisionRes

        indexList = tree.query_ball_point(x=coords, r=radius)

        uniqueIndices = set()
        for indices in indexList:
            uniqueIndices.update(indices)

        n0IndicesCand = numpy.array(list(uniqueIndices))
        trisCand = n0IndicesCand//nNodes
        nextNodeIndex = numpy.mod(n0IndicesCand + 1, nNodes)
        n1IndicesCand = nNodes * trisCand + nextNodeIndex

        n0Cand = Vector(xNode[n0IndicesCand],
                        yNode[n0IndicesCand],
                        zNode[n0IndicesCand])
        n1Cand = Vector(xNode[n1IndicesCand],
                        yNode[n1IndicesCand],
                        zNode[n1IndicesCand])

        intersect = _intersects(n0Cand, n1Cand, transectv0,
                                transectv1)

        n0Inter = Vector(n0Cand.x[intersect],
                         n0Cand.y[intersect],
                         n0Cand.z[intersect])
        n1Inter = Vector(n1Cand.x[intersect],
                         n1Cand.y[intersect],
                         n1Cand.z[intersect])

        trisInter = trisCand[intersect]
        n0IndicesInter = n0IndicesCand[intersect]
        n1IndicesInter = n1IndicesCand[intersect]

        intersections = _intersection(n0Inter, n1Inter, transectv0, transectv1)
        intersections = Vector(earth_radius*intersections.x,
                               earth_radius*intersections.y,
                               earth_radius*intersections.z)

        angularDistance = _angular_distance(first=transectv0,
                                            second=intersections)

        dNodeLocal = dStart + earth_radius * angularDistance

        dStart += earth_radius*_angular_distance(first=transectv0,
                                                 second=transectv1)

        node0Inter = numpy.mod(n0IndicesInter, nNodes)
        node1Inter = numpy.mod(n1IndicesInter, nNodes)

        nodeWeights = (_angular_distance(first=intersections, second=n1Inter) /
                       _angular_distance(first=n0Inter, second=n1Inter))

        weights = numpy.zeros((len(trisInter), nHorizWeights))
        cellIndices = numpy.zeros((len(trisInter), nHorizWeights), int)
        for index in range(3):
            weights[:, index] = (nodeWeights *
                                 nodeCellWeights[trisInter, node0Inter, index])
            cellIndices[:, index] = \
                nodeCellIndices[trisInter, node0Inter, index]
            weights[:, index+3] = ((1.0 - nodeWeights) *
                                   nodeCellWeights[trisInter, node1Inter, index])
            cellIndices[:, index+3] = \
                nodeCellIndices[trisInter, node1Inter, index]

        if first:
            xOut = intersections.x
            yOut = intersections.y
            zOut = intersections.z
            dNode = dNodeLocal

            tris = trisInter
            nodes = node0Inter
            interpCells = cellIndices
            cellWeights = weights
            first = False
        else:
            xOut = numpy.append(xOut, intersections.x)
            yOut = numpy.append(yOut, intersections.y)
            zOut = numpy.append(zOut, intersections.z)
            dNode = numpy.append(dNode, dNodeLocal)

            tris = numpy.concatenate((tris, trisInter))
            nodes = numpy.concatenate((nodes, node0Inter))
            interpCells = numpy.concatenate((interpCells, cellIndices), axis=0)
            cellWeights = numpy.concatenate((cellWeights, weights), axis=0)

        dTransect[segIndex + 1] = dStart

    # sort nodes by distance
    sortIndices = numpy.argsort(dNode)
    dSorted = dNode[sortIndices]
    trisSorted = tris[sortIndices]

    epsilon = 1e-6*subdivisionRes
    nodesAreSame = numpy.abs(dSorted[1:] - dSorted[:-1]) < epsilon
    if nodesAreSame[0]:
        # the first two nodes are the same, so the first transect point is in a
        # triangle, and we need to figure out which
        if trisSorted[1] == trisSorted[2] or trisSorted[1] == trisSorted[3]:
            # the first transect point is in trisSorted[0], so the first two
            # nodes are in the right order
            indices = [0, 1]
        elif trisSorted[0] == trisSorted[2] or trisSorted[0] == trisSorted[3]:
            # the first transect point is in trisSorted[1], so the first two
            # nodes need to be swapped
            indices = [1, 0]
        else:
            raise ValueError('Couldn\'t find an order for the first two nodes')
    else:
        # the first transect point is outside of an MPAS cell
        indices = [0]

    while len(indices) < len(sortIndices):
        index = len(indices)
        currentTri = trisSorted[indices[-1]]
        if nodesAreSame[index]:
            # the next two nodes are the same, so we need to know which
            # corresponds to the current triangle
            if trisSorted[index] == currentTri:
                # the first node is in the current triangle, so add the next
                # two nodes in the current order
                indices.extend([index, index+1])
            elif trisSorted[index+1] == currentTri:
                # the second node is in the current triangle, so add the next
                # two nodes in swapped order
                indices.extend([index+1, index])
            else:
                raise ValueError('Couldn\'t find an order for nodes {} and '
                                 '{}'.format(index, index+1))
        else:
            # the next node is a boundary of the MPAS domain, so there is no
            # ambiguity about order and we just add it
            indices.extend([index])

    indices = sortIndices[indices]

    dNode = dNode[indices]
    xOut= xOut[indices]
    yOut = yOut[indices]
    zOut = zOut[indices]

    tris = tris[indices]
    nodes = nodes[indices]
    interpCells = interpCells[indices, :]
    cellWeights = cellWeights[indices, :]

    # we need to figure out if the end points of the transect are in a triangle
    # and add tris, nodes, interpCells and cellWeights if so

    if len(tris) >= 2 and tris[0] != tris[1]:
        # the starting point is in a triangle so we need to duplicate the first
        # entry in several fields to handle this
        tris = numpy.concatenate((tris[0:1], tris))
        nodes = numpy.concatenate((nodes[0:1], nodes))
        interpCells = numpy.concatenate((interpCells[0:1, :], interpCells),
                                        axis=0)
        cellWeights = numpy.concatenate((cellWeights[0:1, :], cellWeights),
                                        axis=0)

        dNode = numpy.append(numpy.array([0.]), dNode)
        xOut = numpy.append(x[0].values, xOut)
        yOut = numpy.append(y[0].values, yOut)
        zOut = numpy.append(z[0].values, zOut)

    if len(tris) >= 2 and tris[-1] != tris[-2]:
        # the end point is in a triangle so we need to add final entries or
        # duplicate the last entry in several fields to handle this
        tris = numpy.concatenate((tris, tris[-1:]))
        nodes = numpy.concatenate((nodes, nodes[-1:]))
        interpCells = numpy.concatenate((interpCells, interpCells[-1:, :]),
                                        axis=0)
        cellWeights = numpy.concatenate((cellWeights, cellWeights[-1:, :]),
                                        axis=0)

        dNode = numpy.append(dNode, dStart)
        xOut = numpy.append(xOut, x[-1].values)
        yOut = numpy.append(yOut, y[-1].values)
        zOut = numpy.append(zOut, z[-1].values)

    assert(numpy.all(tris[0::2] == tris[1::2]))

    tris = tris[0::2]

    lonOut, latOut = _cartesian_to_lon_lat(xOut, yOut, zOut, earth_radius,
                                           degrees)

    nSegments = len(xOut)//2
    nBounds = 2

    cellIndices = dsTris.triCellIndices.values[tris]
    nodes = nodes.reshape((nSegments, nBounds))
    dNode = dNode.reshape((nSegments, nBounds))

    dsOut = xarray.Dataset()
    dsOut['xCartNode'] = (('nSegments', 'nBounds'),
                          xOut.reshape((nSegments, nBounds)))
    dsOut['yCartNode'] = (('nSegments', 'nBounds'),
                          yOut.reshape((nSegments, nBounds)))
    dsOut['zCartNode'] = (('nSegments', 'nBounds'),
                          zOut.reshape((nSegments, nBounds)))
    dsOut['dNode'] = (('nSegments', 'nBounds'), dNode)
    dsOut['lonNode'] = (('nSegments', 'nBounds'),
                        lonOut.reshape((nSegments, nBounds)))
    dsOut['latNode'] = (('nSegments', 'nBounds'),
                        latOut.reshape((nSegments, nBounds)))

    dsOut['horizTriangleIndices'] = ('nSegments', tris)
    dsOut['cellIndices'] = ('nSegments', cellIndices)
    dsOut['horizTriangleNodeIndices'] = (('nSegments', 'nBounds'), nodes)
    dsOut['interpHorizCellIndices'] = \
        (('nSegments', 'nBounds', 'nHorizWeights'),
         interpCells.reshape((nSegments, nBounds, nHorizWeights)))
    dsOut['interpHorizCellWeights'] = \
        (('nSegments', 'nBounds', 'nHorizWeights'),
         cellWeights.reshape((nSegments, nBounds, nHorizWeights)))

    transectIndicesOnHorizNode = numpy.zeros(dNode.shape, int)
    transectWeightsOnHorizNode = numpy.zeros(dNode.shape)
    for segIndex in range(len(dTransect)-1):
        d0 = dTransect[segIndex]
        d1 = dTransect[segIndex+1]
        mask = numpy.logical_and(dNode >= d0, dNode < d1)
        transectIndicesOnHorizNode[mask] = segIndex
        transectWeightsOnHorizNode[mask] = (d1 - dNode[mask])/(d1 - d0)
    # last index will get missed by the mask and needs to be handled as a
    # special case
    transectIndicesOnHorizNode[-1, 1] = len(dTransect)-2
    transectWeightsOnHorizNode[-1, 1] = 0.0

    dsOut['lonTransect'] = lonTransect
    dsOut['latTransect'] = latTransect
    dsOut['xCartTransect'] = x
    dsOut['yCartTransect'] = y
    dsOut['zCartTransect'] = z
    dsOut['dTransect'] = (lonTransect.dims, dTransect)
    dsOut['transectIndicesOnHorizNode'] = (('nSegments', 'nBounds'),
                                           transectIndicesOnHorizNode)
    dsOut['transectWeightsOnHorizNode'] = (('nSegments', 'nBounds'),
                                           transectWeightsOnHorizNode)

    return dsOut


def subdivide_great_circle(x, y, z, maxRes, earthRadius):  # {{{
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

    angularDistance = _angular_distance(x=x, y=y, z=z)

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

    denom = 1./numpy.sin(delta)
    a = denom*numpy.sin((1.-frac)*delta)
    b = denom*numpy.sin(frac*delta)

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
            earth_radius*_angular_distance(first=transectv0, second=transectv1)

    return distance


def _lon_lat_to_cartesian(lon, lat, earth_radius, degrees):
    """Convert from lon/lat to Cartesian x, y, z"""

    if degrees:
        lon = numpy.deg2rad(lon)
        lat = numpy.deg2rad(lat)
    x = earth_radius * numpy.cos(lat) * numpy.cos(lon)
    y = earth_radius * numpy.cos(lat) * numpy.sin(lon)
    z = earth_radius * numpy.sin(lat)
    return x, y, z


def _cartesian_to_lon_lat(x, y, z, earth_radius, degrees):
    """Convert from  Cartesian x, y, z to lon/lat"""
    lon = numpy.arctan2(y, x)
    lat = numpy.arcsin(z/earth_radius)
    if degrees:
        lon = numpy.rad2deg(lon)
        lat = numpy.rad2deg(lat)
    return lon, lat


def _angular_distance(x=None, y=None, z=None, first=None, second=None):
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


def _intersects(a1, a2, b1, b2):
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


def _intersection(a1, a2, b1, b2):
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
