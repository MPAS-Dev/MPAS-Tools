import xarray
import numpy


def mesh_to_triangles(dsMesh, periodicCopy=False):
    """
    Construct a dataset in which each MPAS cell is divided into the triangles
    connecting pairs of adjacent vertices to cell centers.

    Parameters
    ----------
    dsMesh : xarray.Dataset
        An MPAS mesh

    periodicCopy : bool, optional
        Whether to make a periodic copy of triangles that cross -180/180 degrees
        longitude.  This is helpful when plotting triangles in a lon/lat space.

    Returns
    -------
    dsTris : xarray.Dataset
        A dataset that defines triangles connecting pairs of adjacent vertices
        to cell centers as well as the cell index that each triangle is in and
        cell indices and weights for interpolating data defined at cell centers
        to triangle nodes.  ``dsTris`` includes variables ``triCellIndices``,
        the cell that each triangle is part of; ``nodeCellIndices`` and
        ``nodeCellWeights``, the indices and weights used to interpolate from
        MPAS cell centers to triangle nodes; Cartesian coordinates ``xNode``,
        ``yNode``, and ``zNode``; and ``lonNode``` and ``latNode`` in radians.
        ``lonNode`` is guaranteed to be within 180 degrees of the cell center
        corresponding to ``triCellIndices``.  Nodes always have a
        counterclockwise winding.

    """
    nVerticesOnCell = dsMesh.nEdgesOnCell.values
    verticesOnCell = dsMesh.verticesOnCell.values - 1
    cellsOnVertex = dsMesh.cellsOnVertex.values - 1

    kiteAreasOnVertex = dsMesh.kiteAreasOnVertex.values

    nTriangles = numpy.sum(nVerticesOnCell)

    maxEdges = dsMesh.sizes['maxEdges']
    nCells = dsMesh.sizes['nCells']
    if dsMesh.sizes['vertexDegree'] != 3:
        raise ValueError('mesh_to_triangles only supports meshes with '
                         'vertexDegree = 3')

    # find the third vertex for each triangle
    nextVertex = -1*numpy.ones(verticesOnCell.shape, int)
    for iVertex in range(maxEdges):
        valid = iVertex < nVerticesOnCell
        invalid = numpy.logical_not(valid)
        verticesOnCell[invalid, iVertex] = -1
        nv = nVerticesOnCell[valid]
        cellIndices = numpy.arange(0, nCells)[valid]
        iNext = numpy.where(iVertex < nv - 1, iVertex + 1, 0)
        nextVertex[:, iVertex][valid] = verticesOnCell[cellIndices, iNext]

    valid = verticesOnCell >= 0
    verticesOnCell = verticesOnCell[valid]
    nextVertex = nextVertex[valid]

    # find the cell index for each triangle
    triCellIndices, _ = numpy.meshgrid(numpy.arange(0, nCells),
                                       numpy.arange(0, maxEdges),
                                       indexing='ij')
    triCellIndices = triCellIndices[valid]

    # find list of cells and weights for each triangle node
    nodeCellIndices = -1*numpy.ones((nTriangles, 3, 3), int)
    nodeCellWeights = numpy.zeros((nTriangles, 3, 3))

    # the first node is at the cell center, so the value is just the one from
    # that cell
    nodeCellIndices[:, 0, 0] = triCellIndices
    nodeCellWeights[:, 0, 0] = 1.

    # the other 2 nodes are associated with vertices
    nodeCellIndices[:, 1, :] = cellsOnVertex[verticesOnCell, :]
    nodeCellWeights[:, 1, :] = kiteAreasOnVertex[verticesOnCell, :]
    nodeCellIndices[:, 2, :] = cellsOnVertex[nextVertex, :]
    nodeCellWeights[:, 2, :] = kiteAreasOnVertex[nextVertex, :]

    weightSum = numpy.sum(nodeCellWeights, axis=2)
    for iNode in range(3):
        nodeCellWeights[:, :, iNode] = nodeCellWeights[:, :, iNode]/weightSum

    dsTris = xarray.Dataset()
    dsTris['triCellIndices'] = ('nTriangles', triCellIndices)
    dsTris['nodeCellIndices'] = (('nTriangles', 'nNodes', 'nInterp'),
                                 nodeCellIndices)
    dsTris['nodeCellWeights'] = (('nTriangles', 'nNodes', 'nInterp'),
                                 nodeCellWeights)

    # get Cartesian and lon/lat coordinates of each node
    for prefix in ['x', 'y', 'z', 'lat', 'lon']:
        outVar = '{}Node'.format(prefix)
        cellVar = '{}Cell'.format(prefix)
        vertexVar = '{}Vertex'.format(prefix)
        coord = numpy.zeros((nTriangles, 3))
        coord[:, 0] = dsMesh[cellVar].values[triCellIndices]
        coord[:, 1] = dsMesh[vertexVar].values[verticesOnCell]
        coord[:, 2] = dsMesh[vertexVar].values[nextVertex]
        dsTris[outVar] = (('nTriangles', 'nNodes'), coord)

    # nothing obvious we can do about triangles containing the poles

    # make sure node longitudes are within 180 degrees of the cell center
    lonNode = dsTris.lonNode.values
    lonCell = lonNode[:, 0]
    copyEast = numpy.zeros(lonCell.shape, bool)
    copyWest = numpy.zeros(lonCell.shape, bool)
    for iNode in [1, 2]:
        mask = lonNode[:, iNode] - lonCell > numpy.pi
        copyEast = numpy.logical_or(copyEast, mask)
        lonNode[:, iNode][mask] = lonNode[:, iNode][mask] - 2*numpy.pi
        mask = lonNode[:, iNode] - lonCell < -numpy.pi
        copyWest = numpy.logical_or(copyWest, mask)
        lonNode[:, iNode][mask] = lonNode[:, iNode][mask] + 2*numpy.pi
    if periodicCopy:
        dsNew = xarray.Dataset()
        eastIndices = numpy.nonzero(copyEast)[0]
        westIndices = numpy.nonzero(copyWest)[0]
        triIndices = numpy.append(numpy.append(numpy.arange(0, nTriangles),
                                               eastIndices), westIndices)

        dsNew['triCellIndices'] = ('nTriangles', triCellIndices[triIndices])
        dsNew['nodeCellIndices'] = (('nTriangles', 'nNodes', 'nInterp'),
                                    nodeCellIndices[triIndices, :, :])
        dsNew['nodeCellWeights'] = (('nTriangles', 'nNodes', 'nInterp'),
                                    nodeCellWeights[triIndices, :, :])
        for prefix in ['x', 'y', 'z', 'lat']:
            outVar = '{}Node'.format(prefix)
            coord = dsTris[outVar].values
            dsNew[outVar] = (('nTriangles', 'nNodes'), coord[triIndices, :])

        coord = dsTris['lonNode'].values[triIndices, :]
        eastSlice = slice(nTriangles, nTriangles+len(eastIndices))
        coord[eastSlice, :] = coord[eastSlice, :] + 2*numpy.pi
        westSlice = slice(nTriangles+len(eastIndices),
                          nTriangles + len(eastIndices) + len(westIndices))
        coord[westSlice, :] = coord[westSlice, :] - 2 * numpy.pi
        dsNew['lonNode'] = (('nTriangles', 'nNodes'), coord)
        dsTris = dsNew

    return dsTris
