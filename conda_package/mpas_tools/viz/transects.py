import numpy
import xarray

from scipy.spatial import cKDTree
from shapely.geometry import LineString, Point

from mpas_tools.transects import Vector, lon_lat_to_cartesian, \
    cartesian_to_lon_lat, intersects, intersection, angular_distance, \
    subdivide_great_circle, subdivide_planar


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
        coordinates of these nodes are ``xCartNode``, ``yCartNode``,
        ``zCartNode``, ``lonNode`` and ``latNode``.  The distance along the
        transect of each intersection is ``dNode``. The index of the triangle
        and the first triangle node in ``dsTris`` associated with each
        intersection node are given by ``horizTriangleIndices`` and
        ``horizTriangleNodeIndices``, respectively. The second node on the
        triangle for the edge associated with the intersection is given by
        ``numpy.mod(horizTriangleNodeIndices + 1, 3)``.

        The MPAS cell that a given node belongs to is given by
        ``horizCellIndices``. Each node also has an associated set of 6
        ``interpHorizCellIndices`` and ``interpHorizCellWeights`` that can be
        used to interpolate from MPAS cell centers to nodes first with
        area-weighted averaging to MPAS vertices and then linear interpolation
        along triangle edges.  Some of the weights may be zero, in which case
        the associated ``interpHorizCellIndices`` will be -1.

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

    x, y, z = lon_lat_to_cartesian(lonTransect, latTransect, earth_radius,
                                   degrees)

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

        if len(n0IndicesCand) == 0:
            continue

        trisCand = n0IndicesCand//nNodes
        nextNodeIndex = numpy.mod(n0IndicesCand + 1, nNodes)
        n1IndicesCand = nNodes * trisCand + nextNodeIndex

        n0Cand = Vector(xNode[n0IndicesCand],
                        yNode[n0IndicesCand],
                        zNode[n0IndicesCand])
        n1Cand = Vector(xNode[n1IndicesCand],
                        yNode[n1IndicesCand],
                        zNode[n1IndicesCand])

        intersect = intersects(n0Cand, n1Cand, transectv0,
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

        intersections = intersection(n0Inter, n1Inter, transectv0, transectv1)
        intersections = Vector(earth_radius*intersections.x,
                               earth_radius*intersections.y,
                               earth_radius*intersections.z)

        angularDistance = angular_distance(first=transectv0,
                                            second=intersections)

        dNodeLocal = dStart + earth_radius * angularDistance

        dStart += earth_radius*angular_distance(first=transectv0,
                                                 second=transectv1)

        node0Inter = numpy.mod(n0IndicesInter, nNodes)
        node1Inter = numpy.mod(n1IndicesInter, nNodes)

        nodeWeights = (angular_distance(first=intersections, second=n1Inter) /
                       angular_distance(first=n0Inter, second=n1Inter))

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

    epsilon = 1e-6*subdivisionRes
    dNode, xOut, yOut, zOut, tris, nodes, interpCells, cellWeights = \
        _sort_intersections(dNode, tris, nodes, xOut, yOut, zOut, interpCells,
                            cellWeights, epsilon)

    dNode, xOut, yOut, zOut, tris, nodes, interpCells, cellWeights = \
        _update_start_end_triangles(tris, nodes, interpCells, cellWeights,
                                    dNode, xOut, yOut, zOut, dStart, x, y, z)

    lonOut, latOut = cartesian_to_lon_lat(xOut, yOut, zOut, earth_radius,
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
    dsOut['horizCellIndices'] = ('nSegments', cellIndices)
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


def find_planar_transect_cells_and_weights(xTransect, yTransect, dsTris, dsMesh,
                                           tree, subdivisionRes=10e3):
    """
    Find "nodes" where the transect intersects the edges of the triangles
    that make up MPAS cells.

    Parameters
    ----------
    xTransect, yTransect : xarray.DataArray
        The x and y points defining segments making up the transect

    dsTris : xarray.Dataset
        A dataset that defines triangles, the results of calling
        ``mesh_to_triangles()``

    dsMesh : xarray.Dataset
        A data set with the full MPAS mesh.

    tree : scipy.spatial.cKDTree
        A tree of edge centers from triangles making up an MPAS mesh, the return
        value from ``make_triangle_tree()``

    subdivisionRes : float, optional
        Resolution in m to use to subdivide the transect when looking for
        intersection candidates.  Should be small enough that curvature is
        small.

    Returns
    -------
    dsOut : xarray.Dataset
        A dataset that contains "nodes" where the transect intersects the
        edges of the triangles in ``dsTris``.  The nodes also include the two
        end points of the transect, which typically lie within triangles. Each
        internal node (that is, not including the end points) is purposefully
        repeated twice, once for each triangle that node touches.  This allows
        for discontinuous fields between triangles (e.g. if one wishes to plot
        constant values on each MPAS cell).  The planar coordinates of these
        nodes are ``xNode`` and ``yNode``.  The distance along the transect of
        each intersection is ``dNode``. The index of the triangle and the first
        triangle node in ``dsTris`` associated with each intersection node are
        given by ``horizTriangleIndices`` and ``horizTriangleNodeIndices``,
        respectively. The second node on the triangle for the edge associated
        with the intersection is given by
        ``numpy.mod(horizTriangleNodeIndices + 1, 3)``.

        The MPAS cell that a given node belongs to is given by
        ``horizCellIndices``. Each node also has an associated set of 6
        ``interpHorizCellIndices`` and ``interpHorizCellWeights`` that can be
        used to interpolate from MPAS cell centers to nodes first with
        area-weighted averaging to MPAS vertices and then linear interpolation
        along triangle edges.  Some of the weights may be zero, in which case
        the associated ``interpHorizCellIndices`` will be -1.

        Finally, ``xTransect`` and ``yTransect`` are included in the
        dataset, along with ``dTransect``, the distance along the transect of
        each original transect point.  In order to interpolate values (e.g.
        observations) from the original transect points to the intersection
        nodes, linear interpolation indices ``transectIndicesOnHorizNode`` and
        weights ``transectWeightsOnHorizNode`` are provided.  The values at
        nodes are found by::

          nodeValues = ((transectValues[transectIndicesOnHorizNode] *
                         transectWeightsOnHorizNode)
                        + (transectValues[transectIndicesOnHorizNode+1] *
                           (1.0 - transectWeightsOnHorizNode))
    """
    buffer = numpy.maximum(numpy.amax(dsMesh.dvEdge.values),
                           numpy.amax(dsMesh.dcEdge.values))

    nNodes = dsTris.sizes['nNodes']
    nodeCellWeights = dsTris.nodeCellWeights.values
    nodeCellIndices = dsTris.nodeCellIndices.values

    x = xTransect
    y = yTransect

    xNode = dsTris.xNode.values.ravel()
    yNode = dsTris.yNode.values.ravel()

    coordNode = numpy.zeros((len(xNode), 2))
    coordNode[:, 0] = xNode
    coordNode[:, 1] = yNode

    dTransect = numpy.zeros(xTransect.shape)

    dNode = None
    xOut = None
    yOut = None
    tris = None
    nodes = None
    interpCells = None
    cellWeights = None

    nHorizWeights = 6

    first = True

    dStart = 0.
    for segIndex in range(len(x)-1):

        subSlice = slice(segIndex, segIndex+2)
        xSub, ySub, _, _ = subdivide_planar(
            x[subSlice].values, y[subSlice].values, subdivisionRes)

        startPoint = Point(xTransect[segIndex].values,
                           yTransect[segIndex].values)
        endPoint = Point(xTransect[segIndex+1].values,
                         yTransect[segIndex+1].values)

        segment = LineString([startPoint, endPoint])

        coords = numpy.zeros((len(xSub), 3))
        coords[:, 0] = xSub
        coords[:, 1] = ySub
        radius = buffer + subdivisionRes

        indexList = tree.query_ball_point(x=coords, r=radius)

        uniqueIndices = set()
        for indices in indexList:
            uniqueIndices.update(indices)

        startIndices = numpy.array(list(uniqueIndices))

        if len(startIndices) == 0:
            continue

        trisCand = startIndices//nNodes
        nextNodeIndex = numpy.mod(startIndices + 1, nNodes)
        endIndices = nNodes * trisCand + nextNodeIndex

        intersectingNodes = list()
        trisInter = list()
        xIntersection = list()
        yIntersection = list()
        nodeWeights = list()
        node0Inter = list()
        node1Inter = list()
        distances = list()

        for index in range(len(startIndices)):
            start = startIndices[index]
            end = endIndices[index]

            node0 = Point(coordNode[start, 0], coordNode[start, 1])
            node1 = Point(coordNode[end, 0], coordNode[end, 1])

            edge = LineString([node0, node1])
            if segment.intersects(edge):
                intersection = segment.intersection(edge)
                intersectingNodes.append((node0, node1, start, end, edge))

                if isinstance(intersection, LineString):
                    raise ValueError('A triangle edge exactly coincides with a '
                                     'transect segment and I can\'t handle '
                                     'that case.  Try moving the transect a '
                                     'tiny bit.')
                elif not isinstance(intersection, Point):
                    raise ValueError('Unexpected intersection type {}'.format(
                        intersection))

                xIntersection.append(intersection.x)
                yIntersection.append(intersection.y)

                startToIntersection = LineString([startPoint, intersection])

                weight = (LineString([intersection, node1]).length /
                          LineString([node0, node1]).length)

                nodeWeights.append(weight)
                node0Inter.append(numpy.mod(start, nNodes))
                node1Inter.append(numpy.mod(end, nNodes))
                distances.append(startToIntersection.length)
                trisInter.append(trisCand[index])

        distances = numpy.array(distances)
        xIntersection = numpy.array(xIntersection)
        yIntersection = numpy.array(yIntersection)
        nodeWeights = numpy.array(nodeWeights)
        node0Inter = numpy.array(node0Inter)
        node1Inter = numpy.array(node1Inter)
        trisInter = numpy.array(trisInter)

        dNodeLocal = dStart + distances

        dStart += segment.length

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
            xOut = xIntersection
            yOut = yIntersection
            dNode = dNodeLocal

            tris = trisInter
            nodes = node0Inter
            interpCells = cellIndices
            cellWeights = weights
            first = False
        else:
            xOut = numpy.append(xOut, xIntersection)
            yOut = numpy.append(yOut, yIntersection)
            dNode = numpy.append(dNode, dNodeLocal)

            tris = numpy.concatenate((tris, trisInter))
            nodes = numpy.concatenate((nodes, node0Inter))
            interpCells = numpy.concatenate((interpCells, cellIndices), axis=0)
            cellWeights = numpy.concatenate((cellWeights, weights), axis=0)

        dTransect[segIndex + 1] = dStart

    zOut = numpy.zeros(xOut.shape)
    z = xarray.zeros_like(x)

    epsilon = 1e-6*subdivisionRes
    dNode, xOut, yOut, zOut, tris, nodes, interpCells, cellWeights = \
        _sort_intersections(dNode, tris, nodes, xOut, yOut, zOut, interpCells,
                            cellWeights, epsilon)

    dNode, xOut, yOut, zOut, tris, nodes, interpCells, cellWeights = \
        _update_start_end_triangles(tris, nodes, interpCells, cellWeights,
                                    dNode, xOut, yOut, zOut, dStart, x, y, z)

    nSegments = len(xOut)//2
    nBounds = 2

    cellIndices = dsTris.triCellIndices.values[tris]
    nodes = nodes.reshape((nSegments, nBounds))
    dNode = dNode.reshape((nSegments, nBounds))

    dsOut = xarray.Dataset()
    dsOut['xNode'] = (('nSegments', 'nBounds'),
                      xOut.reshape((nSegments, nBounds)))
    dsOut['yNode'] = (('nSegments', 'nBounds'),
                      yOut.reshape((nSegments, nBounds)))
    dsOut['dNode'] = (('nSegments', 'nBounds'), dNode)

    dsOut['horizTriangleIndices'] = ('nSegments', tris)
    dsOut['horizCellIndices'] = ('nSegments', cellIndices)
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

    dsOut['xTransect'] = x
    dsOut['yTransect'] = y
    dsOut['dTransect'] = (xTransect.dims, dTransect)
    dsOut['transectIndicesOnHorizNode'] = (('nSegments', 'nBounds'),
                                           transectIndicesOnHorizNode)
    dsOut['transectWeightsOnHorizNode'] = (('nSegments', 'nBounds'),
                                           transectWeightsOnHorizNode)

    return dsOut


def _sort_intersections(dNode, tris, nodes, xOut, yOut, zOut, interpCells,
                        cellWeights, epsilon):
    """ sort nodes by distance """

    sortIndices = numpy.argsort(dNode)
    dSorted = dNode[sortIndices]
    trisSorted = tris[sortIndices]

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
        if index < len(nodesAreSame) and nodesAreSame[index]:
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
                print(trisSorted[index:index+2], currentTri)
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

    return dNode, xOut, yOut, zOut, tris, nodes, interpCells, cellWeights


def _update_start_end_triangles(tris, nodes, interpCells, cellWeights, dNode,
                                xOut, yOut, zOut, dStart, x, y, z):
    """
    figure out if the end points of the transect are in a triangle and add tris,
    nodes, interpCells and cellWeights if so
    """

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

    return dNode, xOut, yOut, zOut, tris, nodes, interpCells, cellWeights
