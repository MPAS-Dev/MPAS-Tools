import xarray
import numpy


def find_transect_levels_and_weights(dsTransect, layerThickness, bottomDepth,
                                     maxLevelCell, zTransect=None):
    """
    Construct a vertical coordinate for a transect produced by
    :py:fun:`mpas_tools.viz.transects.find_transect_cells_and_weights()`, then
    break each resulting quadrilateral into 2 triangles that can later be
    visualized with functions like ``tripcolor`` and ``tricontourf``.  Also,
    compute interpolation weights such that observations at points on the
    original transect and with vertical coordinate ``transectZ`` can be
    bilinearly interpolated to the nodes of the transect.

    Parameters
    ----------
    dsTransect : xarray.Dataset
        A dataset that defines nodes of the transect, the results of calling
        :py:fun:`mpas_tools.viz.transects.find_transect_cells_and_weights()`

    layerThickness : xarray.DataArray
        layer thicknesses on the MPAS mesh

    bottomDepth : xarray.DataArray
        the (positive down) depth of the seafloor on the MPAS mesh

    maxLevelCell : xarray.DataArray
        the vertical zero-based index of the bathymetry on the MPAS mesh

    zTransect : xarray.DataArray, optional
        the z coordinate of the transect (1D or 2D).  If 2D, it must have the
        same along-transect dimension as the lon and lat passed to
        :py:fun:`mpas_tools.viz.transects.find_transect_cells_and_weights()`

    Returns
    -------
    dsTransectTriangles : xarray.Dataset
        A dataset that contains nodes and triangles that make up a 2D transect.
        For convenience in visualization, the quadrilaterals of each cell making
        up the transect have been divided into an upper and lower triangle.  The
        nodes of the triangles are completely independent of one another,
        allowing for potential jumps in fields values between nodes of different
        triangles that are at the same location.  This is convenient, for
        example, when visualizing data with constant values within each MPAS
        cell.

        There are ``nTransectTriangles = 2*nTransectCells`` triangles, each with
        ``nTriangleNodes = 3`` nodes, where ``nTransectCells`` is the number of
        valid transect cells (quadrilaterals) that are above the MPAS-Ocean
        bathymetry.

        In addition to the variables and coordinates in ``dsTransect``, the
        output dataset contains:

            - nodeHorizBoundsIndices: which of the ``nHorizBounds = 2``
              bounds of a horizontal transect segment a given node is on
            - segmentIndices: the transect segment of a given triangle
            - cellIndices: the MPAS-Ocean cell of a given triangle
            - levelIndices: the MPAS-Ocean vertical level of a given triangle

            - zTransectNode: the vertical height of each triangle node
            - ssh, zSeaFloor: the sea-surface height and sea-floor height at
              each node of each transect segment

            - interpCellIndices, interpLevelIndices: the MPAS-Ocean cells and
              levels from which the value at a given triangle node are
              interpolated.  This can involve up to ``nWeights = 12`` different
              cells and levels.
            - interpCellWeights: the weight to multiply each field value by
              to perform interpolation to a triangle node.

            - transectInterpVertIndices, transectInterpVertWeights - if
              ``zTransect`` is not ``None``, vertical indices and weights for
              interpolating from the original transect grid to MPAS-Ocean
              transect nodes.

        Interpolation of a DataArray from MPAS cells and levels to transect
        triangle nodes can be performed with
        ``interp_mpas_to_transect_triangle_nodes()``.  Similarly, interpolation of a
        DataArray (e.g. an observation) from the original transect grid to
        transect triangle nodes can be performed with
        ``interp_transect_grid_to_transect_triangle_nodes()``

        To visualize constant values on MPAS cells, a field can be sampled
        at indices ``nCells=cellIndices`` and ``nVertLevels=levelIndices``.
        If a smoother visualization is desired, bilinear interpolation can be
        performed by first sampling the field at indices
        ``nCells=interpCellIndices`` and ``nVertLevels=interpLevelIndices`` and
        then multiplying by ``interpCellWeights`` and summing over
        ``nWeights``.
    """

    dsTransectCells = _get_transect_cells_and_nodes(
        dsTransect, layerThickness, bottomDepth, maxLevelCell)

    dsTransectTriangles = _transect_cells_to_triangles(dsTransectCells)

    if zTransect is not None:
        dsTransectTriangles = _add_vertical_interpolation_of_transect_points(
            dsTransectTriangles, zTransect)

    return dsTransectTriangles


def interp_mpas_to_transect_triangles(dsTransectTriangles, da):
    """
    Interpolate a 3D (``nCells`` by ``nVertLevels``) MPAS-Ocean DataArray
    to transect nodes with constant values in each MPAS cell

    Parameters
    ----------
    dsTransectTriangles : xarray.Dataset
        A dataset that defines triangles making up an MPAS-Ocean transect, the
        results of calling ``find_transect_levels_and_weights()``

    da : xarray.DataArray
        An MPAS-Ocean 3D field with dimensions `nCells`` and ``nVertLevels``
        (possibly among others)

    Returns
    -------
    daNodes : xarray.DataArray
        The data array interpolated to transect nodes with dimensions
        ``nTransectTriangles`` and ``nTriangleNodes`` (in addition to whatever
        dimensions were in ``da`` besides ``nCells`` and ``nVertLevels``)
    """

    cellIndices = dsTransectTriangles.cellIndices
    levelIndices = dsTransectTriangles.levelIndices

    daNodes = da.isel(nCells=cellIndices, nVertLevels=levelIndices)

    return daNodes


def interp_mpas_to_transect_triangle_nodes(dsTransectTriangles, da):
    """
    Interpolate a 3D (``nCells`` by ``nVertLevels``) MPAS-Ocean DataArray
    to transect nodes, linearly interpolating fields between the closest
    neighboring cells

    Parameters
    ----------
    dsTransectTriangles : xarray.Dataset
        A dataset that defines triangles making up an MPAS-Ocean transect, the
        results of calling ``find_transect_levels_and_weights()``

    da : xarray.DataArray
        An MPAS-Ocean 3D field with dimensions `nCells`` and ``nVertLevels``
        (possibly among others)

    Returns
    -------
    daNodes : xarray.DataArray
        The data array interpolated to transect nodes with dimensions
        ``nTransectTriangles`` and ``nTriangleNodes`` (in addition to whatever
        dimensions were in ``da`` besides ``nCells`` and ``nVertLevels``)
    """

    interpCellIndices = dsTransectTriangles.interpCellIndices
    interpLevelIndices = dsTransectTriangles.interpLevelIndices
    interpCellWeights = dsTransectTriangles.interpCellWeights

    da = da.isel(nCells=interpCellIndices, nVertLevels=interpLevelIndices)

    daNodes = (da*interpCellWeights).sum(dim='nWeights')

    return daNodes


def interp_transect_grid_to_transect_triangle_nodes(dsTransectTriangles, da):
    """
    Interpolate a DataArray on the original transect grid to triangle nodes on
    the MPAS-Ocean transect.

    Parameters
    ----------
    dsTransectTriangles : xarray.Dataset
        A dataset that defines triangles making up an MPAS-Ocean transect, the
        results of calling ``find_transect_levels_and_weights()``

    da : xarray.DataArray
        An field on the original triangle mesh

    Returns
    -------
    daNodes : xarray.DataArray
        The data array interpolated to transect nodes with dimensions
        ``nTransectTriangles`` and ``nTriangleNodes``
    """

    horizDim = dsTransectTriangles.dTransect.dims[0]
    zTransect = dsTransectTriangles.zTransect
    vertDim = None
    for dim in zTransect.dims:
        if dim != horizDim:
            vertDim = dim

    horizIndices = dsTransectTriangles.transectIndicesOnHorizNode
    horizWeights = dsTransectTriangles.transectWeightsOnHorizNode

    segmentIndices = dsTransectTriangles.segmentIndices
    nodeHorizBoundsIndices = dsTransectTriangles.nodeHorizBoundsIndices

    horizIndices = horizIndices.isel(nSegments=segmentIndices,
                                     nHorizBounds=nodeHorizBoundsIndices)
    horizWeights = horizWeights.isel(nSegments=segmentIndices,
                                     nHorizBounds=nodeHorizBoundsIndices)

    vertIndices = dsTransectTriangles.transectInterpVertIndices
    vertWeights = dsTransectTriangles.transectInterpVertWeights

    kwargs00 = {horizDim: horizIndices, vertDim: vertIndices}
    kwargs01 = {horizDim: horizIndices, vertDim: vertIndices+1}
    kwargs10 = {horizDim: horizIndices+1, vertDim: vertIndices}
    kwargs11 = {horizDim: horizIndices+1, vertDim: vertIndices+1}

    daNodes = (horizWeights * vertWeights * da.isel(**kwargs00) +
               horizWeights * (1.0 - vertWeights) * da.isel(**kwargs01) +
               (1.0 - horizWeights) * vertWeights * da.isel(**kwargs10) +
               (1.0 - horizWeights) * (1.0 - vertWeights) * da.isel(**kwargs11))

    mask = numpy.logical_and(horizIndices != -1, vertIndices != -1)

    daNodes = daNodes.where(mask)

    return daNodes


def get_outline_segments(dsTransectTriangles, epsilon=1e-3):
    """Get a set of line segments that outline the given transect"""

    nSegments = dsTransectTriangles.sizes['nSegments']
    dSeaFloor = dsTransectTriangles.dNode.values
    zSeaFloor = dsTransectTriangles.zSeaFloor.values
    ssh = dsTransectTriangles.ssh.values

    seaFloorJumps = numpy.abs(dSeaFloor[0:-1, 1] - dSeaFloor[1:, 0]) < epsilon
    nSeaFloorJumps = numpy.count_nonzero(seaFloorJumps)
    nSurface = nSegments
    nSeaFloor = nSegments + nSeaFloorJumps
    nLandJumps = nSegments - nSeaFloorJumps

    nOutline = nSurface + nSeaFloor + 2*nLandJumps

    d = numpy.zeros((nOutline, 2))
    z = numpy.zeros((nOutline, 2))

    d[0:nSegments, :] = dSeaFloor
    z[0:nSegments, :] = ssh
    d[nSegments:2*nSegments, :] = dSeaFloor
    z[nSegments:2*nSegments, :] = zSeaFloor

    dSeaFloorJump = numpy.vstack((dSeaFloor[0:-1, 1], dSeaFloor[1:, 0])).T
    zSeaFloorJump = numpy.vstack((zSeaFloor[0:-1, 1], zSeaFloor[1:, 0])).T
    slc = slice(2*nSegments, 2*nSegments+nSeaFloorJumps)
    d[slc, :] = dSeaFloorJump[seaFloorJumps, :]
    z[slc, :] = zSeaFloorJump[seaFloorJumps, :]

    landJumps1 = numpy.ones(nSegments, bool)
    landJumps1[1:] = numpy.logical_not(seaFloorJumps)
    landJumps2 = numpy.ones(nSegments, bool)
    landJumps2[0:-1] = numpy.logical_not(seaFloorJumps)

    offset = 2*nSegments+nSeaFloorJumps
    slc = slice(offset, offset + nLandJumps)
    d[slc, 0] = dSeaFloor[landJumps1, 0]
    d[slc, 1] = dSeaFloor[landJumps1, 0]
    z[slc, 0] = ssh[landJumps1, 0]
    z[slc, 1] = zSeaFloor[landJumps1, 0]

    offset = 2*nSegments+nSeaFloorJumps+nLandJumps
    slc = slice(offset, offset + nLandJumps)
    d[slc, 0] = dSeaFloor[landJumps2, 1]
    d[slc, 1] = dSeaFloor[landJumps2, 1]
    z[slc, 0] = ssh[landJumps2, 1]
    z[slc, 1] = zSeaFloor[landJumps2, 1]

    d = d.T
    z = z.T

    return d, z


def _get_transect_cells_and_nodes(dsTransect, layerThickness, bottomDepth,
                                  maxLevelCell):

    if 'Time' in layerThickness.dims:
        raise ValueError('Please select a single time level in layerThickness.')

    dsTransect = dsTransect.rename({'nBounds': 'nHorizBounds'})

    zTop, zMid, zBot, ssh, zSeaFloor, interpCellIndices, interpCellWeights = \
        _get_vertical_coordinate(dsTransect, layerThickness, bottomDepth,
                                 maxLevelCell)

    nVertLevels = layerThickness.sizes['nVertLevels']

    levelIndices = xarray.DataArray(data=numpy.arange(2*nVertLevels)//2,
                                    dims='nHalfLevels')

    cellMask = (levelIndices <= maxLevelCell).transpose('nCells', 'nHalfLevels')

    dsTransectCells = _add_valid_cells_and_levels(
        dsTransect, dsTransect.horizCellIndices.values, levelIndices.values,
        cellMask.values)

    # transect cells are made up of half-levels, and each half-level has a top
    # and bottom interface, so we need 4 interfaces per MPAS level

    interpCellIndices, interpLevelIndices, interpCellWeights = \
        _get_interp_indices_and_weights(layerThickness, maxLevelCell,
                                        interpCellIndices, interpCellWeights)

    levelIndex, tempIndex = numpy.meshgrid(numpy.arange(nVertLevels),
                                           numpy.arange(2), indexing='ij')
    levelIndex = xarray.DataArray(data=levelIndex.ravel(), dims='nHalfLevels')
    tempIndex = xarray.DataArray(data=tempIndex.ravel(), dims='nHalfLevels')
    zTop = xarray.concat((zTop, zMid), dim='nTemp')
    zTop = zTop.isel(nVertLevels=levelIndex, nTemp=tempIndex)
    zBot = xarray.concat((zMid, zBot), dim='nTemp')
    zBot = zBot.isel(nVertLevels=levelIndex, nTemp=tempIndex)

    zInterface = xarray.concat((zTop, zBot), dim='nVertBounds')

    segmentIndices = dsTransectCells.segmentIndices
    halfLevelIndices = dsTransectCells.halfLevelIndices

    dsTransectCells['interpCellIndices'] = interpCellIndices.isel(
        nSegments=segmentIndices, nHalfLevels=halfLevelIndices)
    dsTransectCells['interpLevelIndices'] = interpLevelIndices.isel(
        nSegments=segmentIndices, nHalfLevels=halfLevelIndices)
    dsTransectCells['interpCellWeights'] = interpCellWeights.isel(
        nSegments=segmentIndices, nHalfLevels=halfLevelIndices)
    dsTransectCells['zTransectNode'] = zInterface.isel(
        nSegments=segmentIndices, nHalfLevels=halfLevelIndices)

    dsTransectCells['ssh'] = ssh
    dsTransectCells['zSeaFloor'] = zSeaFloor

    dims = ['nSegments', 'nTransectCells', 'nHorizBounds', 'nVertBounds',
            'nHorizWeights', 'nWeights']
    for dim in dsTransectCells.dims:
        if dim not in dims:
            dims.insert(0, dim)
    dsTransectCells = dsTransectCells.transpose(*dims)

    return dsTransectCells


def _get_vertical_coordinate(dsTransect, layerThickness, bottomDepth,
                             maxLevelCell):
    nVertLevels = layerThickness.sizes['nVertLevels']
    levelIndices = xarray.DataArray(data=numpy.arange(nVertLevels),
                                    dims='nVertLevels')
    cellMask = (levelIndices <= maxLevelCell).transpose('nCells', 'nVertLevels')

    ssh = -bottomDepth + layerThickness.sum(dim='nVertLevels')

    interpCellIndices = dsTransect.interpHorizCellIndices
    interpCellWeights = dsTransect.interpHorizCellWeights

    interpMask = numpy.logical_and(interpCellIndices > 0,
                                   cellMask.isel(nCells=interpCellIndices))

    interpCellWeights = interpMask*interpCellWeights
    weightSum = interpCellWeights.sum(dim='nHorizWeights')

    cellIndices = dsTransect.horizCellIndices

    validCells = cellMask.isel(nCells=cellIndices)

    _, validWeights = xarray.broadcast(interpCellWeights, validCells)
    interpCellWeights = (interpCellWeights/weightSum).where(validWeights)

    layerThicknessTransect = layerThickness.isel(nCells=interpCellIndices)
    layerThicknessTransect = (layerThicknessTransect*interpCellWeights).sum(
        dim='nHorizWeights')

    sshTransect = ssh.isel(nCells=interpCellIndices)
    sshTransect = (sshTransect*dsTransect.interpHorizCellWeights).sum(
        dim='nHorizWeights')

    zBot = sshTransect - layerThicknessTransect.cumsum(dim='nVertLevels')
    zTop = zBot + layerThicknessTransect
    zMid = 0.5*(zTop + zBot)

    zSeaFloor = sshTransect - layerThicknessTransect.sum(dim='nVertLevels')

    return zTop, zMid, zBot, sshTransect, zSeaFloor, interpCellIndices, \
        interpCellWeights


def _add_valid_cells_and_levels(dsTransect, cellIndices, levelIndices,
                                cellMask):

    dims = ('nTransectCells',)
    CellIndices, LevelIndices = numpy.meshgrid(cellIndices, levelIndices,
                                               indexing='ij')
    mask = numpy.logical_and(CellIndices >= 0, cellMask[cellIndices, :])

    SegmentIndices, HalfLevelIndices = \
        numpy.meshgrid(numpy.arange(len(cellIndices)),
                       numpy.arange(len(levelIndices)), indexing='ij')

    segmentIndices = xarray.DataArray(data=SegmentIndices[mask], dims=dims)

    dsTransectCells = dsTransect
    dsTransectCells['cellIndices'] = (dims, CellIndices[mask])
    dsTransectCells['levelIndices'] = (dims, LevelIndices[mask])
    dsTransectCells['segmentIndices'] = segmentIndices
    dsTransectCells['halfLevelIndices'] = (dims, HalfLevelIndices[mask])

    return dsTransectCells


def _get_interp_indices_and_weights(layerThickness, maxLevelCell,
                                    interpCellIndices, interpCellWeights):
    interpCellIndices = interpCellIndices.rename({'nHorizWeights': 'nWeights'})
    interpCellWeights = interpCellWeights.rename({'nHorizWeights': 'nWeights'})
    nVertLevels = layerThickness.sizes['nVertLevels']
    nHalfLevels = 2*nVertLevels
    nVertBounds = 2

    interpMaxLevelCell = maxLevelCell.isel(nCells=interpCellIndices)

    levelIndices = xarray.DataArray(
        data=numpy.arange(nHalfLevels)//2, dims='nHalfLevels')
    valid = levelIndices <= interpMaxLevelCell

    topLevelIndices = -1*numpy.ones((nHalfLevels, nVertBounds), int)
    topLevelIndices[1:, 0] = numpy.arange(nHalfLevels-1)//2
    topLevelIndices[:, 1] = numpy.arange(nHalfLevels)//2
    topLevelIndices = xarray.DataArray(
        data=topLevelIndices, dims=('nHalfLevels', 'nVertBounds'))
    interpCellIndices, topLevelIndices = \
        xarray.broadcast(interpCellIndices, topLevelIndices)
    topLevelIndices = topLevelIndices.where(valid, -1)

    botLevelIndices = numpy.zeros((nHalfLevels, nVertBounds), int)
    botLevelIndices[:, 0] = numpy.arange(nHalfLevels)//2
    botLevelIndices[:, 1] = numpy.arange(1, nHalfLevels+1)//2
    botLevelIndices = xarray.DataArray(
        data=botLevelIndices, dims=('nHalfLevels', 'nVertBounds'))
    _, botLevelIndices = xarray.broadcast(interpCellIndices, botLevelIndices)
    botLevelIndices = botLevelIndices.where(valid, -1)
    botLevelIndices = numpy.minimum(botLevelIndices, interpMaxLevelCell)

    interpLevelIndices = xarray.concat((topLevelIndices, botLevelIndices),
                                       dim='nWeights')

    topHalfLevelThickness = 0.5*layerThickness.isel(
        nCells=interpCellIndices, nVertLevels=topLevelIndices)
    topHalfLevelThickness = topHalfLevelThickness.where(topLevelIndices >= 0,
                                                        other=0.)
    botHalfLevelThickness = 0.5*layerThickness.isel(
        nCells=interpCellIndices, nVertLevels=botLevelIndices)

    # vertical weights are proportional to the half-level thickness
    interpCellWeights = xarray.concat(
        (topHalfLevelThickness*interpCellWeights.isel(nVertLevels=topLevelIndices),
         botHalfLevelThickness*interpCellWeights.isel(nVertLevels=botLevelIndices)),
        dim='nWeights')

    weightSum = interpCellWeights.sum(dim='nWeights')
    _, outMask = xarray.broadcast(interpCellWeights, weightSum > 0.)
    interpCellWeights = (interpCellWeights/weightSum).where(outMask)

    interpCellIndices = xarray.concat((interpCellIndices, interpCellIndices),
                                      dim='nWeights')

    return interpCellIndices, interpLevelIndices, interpCellWeights


def _transect_cells_to_triangles(dsTransectCells):

    nTransectCells = dsTransectCells.sizes['nTransectCells']
    nTransectTriangles = 2*nTransectCells
    triangleTransectCellIndices = numpy.arange(nTransectTriangles)//2
    nodeTransectCellIndices = numpy.zeros((nTransectTriangles, 3), int)
    nodeHorizBoundsIndices = numpy.zeros((nTransectTriangles, 3), int)
    nodeVertBoundsIndices = numpy.zeros((nTransectTriangles, 3), int)

    for index in range(3):
        nodeTransectCellIndices[:, index] = triangleTransectCellIndices

    # the upper triangle
    nodeHorizBoundsIndices[0::2, 0] = 0
    nodeVertBoundsIndices[0::2, 0] = 0
    nodeHorizBoundsIndices[0::2, 1] = 1
    nodeVertBoundsIndices[0::2, 1] = 0
    nodeHorizBoundsIndices[0::2, 2] = 0
    nodeVertBoundsIndices[0::2, 2] = 1

    # the lower triangle
    nodeHorizBoundsIndices[1::2, 0] = 0
    nodeVertBoundsIndices[1::2, 0] = 1
    nodeHorizBoundsIndices[1::2, 1] = 1
    nodeVertBoundsIndices[1::2, 1] = 0
    nodeHorizBoundsIndices[1::2, 2] = 1
    nodeVertBoundsIndices[1::2, 2] = 1

    triangleTransectCellIndices = xarray.DataArray(
        data=triangleTransectCellIndices, dims='nTransectTriangles')
    nodeTransectCellIndices = xarray.DataArray(
        data=nodeTransectCellIndices,
        dims=('nTransectTriangles', 'nTriangleNodes'))
    nodeHorizBoundsIndices = xarray.DataArray(
        data=nodeHorizBoundsIndices,
        dims=('nTransectTriangles', 'nTriangleNodes'))
    nodeVertBoundsIndices = xarray.DataArray(
        data=nodeVertBoundsIndices,
        dims=('nTransectTriangles', 'nTriangleNodes'))

    dsTransectTriangles = xarray.Dataset()
    dsTransectTriangles['nodeHorizBoundsIndices'] = \
        nodeHorizBoundsIndices
    for var_name in dsTransectCells.data_vars:
        var = dsTransectCells[var_name]
        if 'nTransectCells' in var.dims:
            if 'nVertBounds' in var.dims:
                assert 'nHorizBounds' in var.dims
                dsTransectTriangles[var_name] = var.isel(
                    nTransectCells=nodeTransectCellIndices,
                    nHorizBounds=nodeHorizBoundsIndices,
                    nVertBounds=nodeVertBoundsIndices)
            elif 'nHorizBounds' in var.dims:
                dsTransectTriangles[var_name] = var.isel(
                    nTransectCells=nodeTransectCellIndices,
                    nHorizBounds=nodeHorizBoundsIndices)
            else:
                dsTransectTriangles[var_name] = var.isel(
                    nTransectCells=triangleTransectCellIndices)
        else:
            dsTransectTriangles[var_name] = var

    dsTransectTriangles = dsTransectTriangles.drop_vars('halfLevelIndices')

    return dsTransectTriangles


def _add_vertical_interpolation_of_transect_points(dsTransectTriangles,
                                                   zTransect):

    dTransect = dsTransectTriangles.dTransect
    # make sure zTransect is 2D
    zTransect, _ = xarray.broadcast(zTransect, dTransect)

    assert len(zTransect.dims) == 2

    horizDim = dTransect.dims[0]
    vertDim = None
    for dim in zTransect.dims:
        if dim != horizDim:
            vertDim = dim

    assert vertDim is not None

    nzTransect = zTransect.sizes[vertDim]

    horizIndices = dsTransectTriangles.transectIndicesOnHorizNode
    horizWeights = dsTransectTriangles.transectWeightsOnHorizNode
    kwargs0 = {horizDim: horizIndices}
    kwargs1 = {horizDim: horizIndices+1}
    zTransectAtHorizNodes = horizWeights*zTransect.isel(**kwargs0) + \
                            (1.0 - horizWeights)*zTransect.isel(**kwargs1)

    zTriangleNode = dsTransectTriangles.zTransectNode

    segmentIndices = dsTransectTriangles.segmentIndices
    nodeHorizBoundsIndices = dsTransectTriangles.nodeHorizBoundsIndices

    nTransectTriangles = dsTransectTriangles.sizes['nTransectTriangles']
    nTriangleNodes = dsTransectTriangles.sizes['nTriangleNodes']
    transectInterpVertIndices = -1*numpy.ones(
        (nTransectTriangles, nTriangleNodes), int)
    transectInterpVertWeights = numpy.zeros(
        (nTransectTriangles, nTriangleNodes))

    kwargs = {vertDim: 0, 'nSegments': segmentIndices,
              'nHorizBounds': nodeHorizBoundsIndices}
    z0 = zTransectAtHorizNodes.isel(**kwargs)
    for zIndex in range(nzTransect-1):
        kwargs = {vertDim: zIndex+1, 'nSegments': segmentIndices,
                  'nHorizBounds': nodeHorizBoundsIndices}
        z1 = zTransectAtHorizNodes.isel(**kwargs)
        mask = numpy.logical_and(zTriangleNode <= z0, zTriangleNode > z1)
        mask = mask.values
        weights = (z1 - zTriangleNode)/(z1 - z0)

        transectInterpVertIndices[mask] = zIndex
        transectInterpVertWeights[mask] = weights.values[mask]
        z0 = z1

    dsTransectTriangles['transectInterpVertIndices'] = (
        ('nTransectTriangles', 'nTriangleNodes'), transectInterpVertIndices)

    dsTransectTriangles['transectInterpVertWeights'] = (
        ('nTransectTriangles', 'nTriangleNodes'), transectInterpVertWeights)

    dsTransectTriangles['zTransect'] = zTransect

    return dsTransectTriangles
