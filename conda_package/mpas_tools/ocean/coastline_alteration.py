from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import xarray


def add_critical_land_blockages(dsMask, dsBlockages):
    '''
    Parameters
    ----------
    dsMask : `xarray.Dataset`
        The mask to which critical blockages should be added
    dsBlockage : `xarray.Dataset`
        The transect masks defining critical land regions that should block
        ocean flow (e.g. the Antarctic Peninsula)

    Returns
    -------
    dsMask : `xarray.Dataset`
        The mask with critical blockages included
    '''

    dsMask = dsMask.copy()

    nTransects = dsBlockages.sizes['nTransects']
    for transectIndex in range(nTransects):
        dsMask.regionCellMasks[:, 0] = numpy.maximum(
            dsBlockages.transectCellMasks[:, transectIndex],
            dsMask.regionCellMasks[:, 0])

    return dsMask


def widen_transect_edge_masks(dsMask, dsMesh, latitude_threshold=43.0):
    '''
    Parameters
    ----------
    dsMask : `xarray.Dataset`
        The mask to which critical blockages should be added
    dsMesh : `xarray.Dataset`
        The transect masks defining critical land regions that should block
        ocean flow (e.g. the Antarctic Peninsula)
    latitude_threshold : float
        Minimum latitude, degrees, for transect widening

    Returns
    -------
    dsMask : `xarray.Dataset`
        The mask with critical blockages included
    '''
    latitude_threshold_radians = numpy.deg2rad(latitude_threshold)

    dsMask = dsMask.copy()

    maxEdges = dsMesh.sizes['maxEdges']

    latMask = numpy.abs(dsMesh.latEdge) > latitude_threshold_radians

    edgeMask = numpy.logical_and(
        latMask, dsMask.transectEdgeMasks == 1)
    for iEdge in range(maxEdges):
        eoc = dsMesh.edgesOnCell[:, iEdge]-1
        mask = numpy.logical_and(eoc >= 0,
                                 edgeMask[eoc])
        # cells with a neighboring transect edge should be masked to 1
        dsMask['transectCellMasks'] = dsMask.transectCellMasks.where(
            numpy.logical_not(mask), 1.)

    return dsMask


def add_land_locked_cells_to_mask(dsMask, dsMesh, latitude_threshold=43.0,
                                  nSweeps=10):
    '''
    Find ocean cells that are land-locked, and alter the cell mask so that they
    are counted as land cells.

    Parameters
    ----------
    dsMask : ``xarray.Dataset``
        A land-mask data set

    dsMesh : ``xarray.Dataset``
        MPAS Mesh data set

    latitude_threshold : float, optional
        Minimum latitude, in degrees, for transect widening

    nSweeps : int, optional
        Maximum number of sweeps to search for land-locked cells

    Returns
    -------
    dsMask : ``xarray.Dataset``
        A copy of the land-mask data set with land-locked cells added to the
        mask for the first region
    '''

    dsMask = xarray.Dataset(dsMask)
    dsMesh = dsMesh.copy(deep=True)

    landMask = dsMask.regionCellMasks.max(dim='nRegions') > 0

    dsMask['landMaskDiagnostic'] = xarray.where(landMask, 1, 0)

    print("Running add_land_locked_cells_to_mask.py.  Total number of cells: "
          "{}".format(dsMesh.sizes['nCells']))

    cellsOnCell = dsMesh.cellsOnCell - 1
    nEdgesOnCell = dsMesh.nEdgesOnCell

    nextCellsOnCell = cellsOnCell.copy(deep=True)
    prevCellsOnCell = cellsOnCell.copy(deep=True)
    for iEdgeOnCell in range(nextCellsOnCell.shape[1]):
        iP1 = numpy.mod(iEdgeOnCell + 1, nEdgesOnCell)
        nextCellsOnCell[:, iEdgeOnCell] = cellsOnCell[:, iP1]
        iM1 = numpy.mod(iEdgeOnCell - 1, nEdgesOnCell)
        prevCellsOnCell[:, iEdgeOnCell] = cellsOnCell[:, iM1]

    dsMesh['cellsOnCell'] = cellsOnCell
    dsMesh['nextCellsOnCell'] = nextCellsOnCell
    dsMesh['prevCellsOnCell'] = prevCellsOnCell
    dsMesh['latCell'] = numpy.rad2deg(dsMesh.latCell)
    dsMesh['lonCell'] = numpy.rad2deg(dsMesh.lonCell)

    landMask, removable = _remove_cells_with_isolated_edges1(
        dsMask, dsMesh, landMask, latitude_threshold)
    landMask = _remove_cells_with_isolated_edges2(
        dsMask, dsMesh, landMask, removable, nSweeps)
    oceanMask = _flood_fill(dsMask, dsMesh, landMask, removable)
    landMask = _revert_cells_with_connected_edges(
        dsMask, dsMesh, oceanMask, landMask, removable, nSweeps)

    return dsMask


def _remove_cells_with_isolated_edges1(dsMask, dsMesh, landMask,
                                       latitude_threshold):
    print("Step 1: Searching for land-locked cells.  Remove cells that only "
          "have isolated active edges.")

    landMaskNew = landMask.copy(deep=True)

    active = numpy.logical_not(landMask)
    removable = numpy.logical_and(
        numpy.abs(dsMesh.latCell) >= latitude_threshold, active)

    cellsOnCell = dsMesh.cellsOnCell
    valid = numpy.logical_and(removable, cellsOnCell >= 0)
    activeEdge = numpy.logical_and(valid, active[cellsOnCell])

    nextCellsOnCell = dsMesh.nextCellsOnCell
    valid = numpy.logical_and(removable, nextCellsOnCell >= 0)
    activeNextEdge = numpy.logical_and(valid, active[nextCellsOnCell])

    # which vertices have adjacent active edges on this cell?
    activeAdjacentEdges = numpy.logical_and(activeEdge, activeNextEdge)

    # which removable cells have no pairs of adjacent active cells?
    noActiveAdjacentEdges = numpy.logical_and(
        removable, numpy.logical_not(numpy.any(activeAdjacentEdges, axis=1)))

    landMaskNew[noActiveAdjacentEdges] = 1
    landLockedCounter = numpy.count_nonzero(noActiveAdjacentEdges)

    dsMask.regionCellMasks[:, 0] = numpy.maximum(dsMask.regionCellMasks[:, 0],
                                                 1*noActiveAdjacentEdges)

    dsMask.landMaskDiagnostic[noActiveAdjacentEdges] = 2

    print("  Number of landLocked cells: {}".format(landLockedCounter))

    return landMaskNew, removable


def _remove_cells_with_isolated_edges2(dsMask, dsMesh, landMask, removable,
                                       nSweeps):
    print("Step 2: Searching for land-locked cells. Remove cells that have "
          "any isolated active edges.")

    cellsOnCell = dsMesh.cellsOnCell
    nextCellsOnCell = dsMesh.nextCellsOnCell
    prevCellsOnCell = dsMesh.prevCellsOnCell

    for iSweep in range(nSweeps):
        landLockedCounter = 0
        landMaskNew = landMask.copy(deep=True)

        active = numpy.logical_not(landMask)
        mask = numpy.logical_and(removable, active)

        valid = numpy.logical_and(mask, cellsOnCell >= 0)
        activeEdge = numpy.logical_and(valid, active[cellsOnCell])
        valid = numpy.logical_and(mask, nextCellsOnCell >= 0)
        activeNextEdge = numpy.logical_and(valid, active[nextCellsOnCell])
        valid = numpy.logical_and(mask, prevCellsOnCell >= 0)
        activePrevEdge = numpy.logical_and(valid, active[prevCellsOnCell])

        # an edge is land-locked if it is active but neither neighbor is active
        landLockedEdges = numpy.logical_and(
            activeEdge,
            numpy.logical_not(
                numpy.logical_or(activePrevEdge, activeNextEdge)))

        landLockedCells = numpy.any(landLockedEdges, axis=1)

        landLockedCounter = numpy.count_nonzero(landLockedCells)
        if landLockedCounter > 0:
            landMaskNew[landLockedCells] = 1
            dsMask.regionCellMasks[landLockedCells, 0] = 1
            dsMask.landMaskDiagnostic[landLockedCells] = 3

        landMask = landMaskNew
        print("  Sweep: {} Number of landLocked cells removed: {}".format(
            iSweep + 1, landLockedCounter))
        if landLockedCounter == 0:
            break

    return landMask


def _flood_fill(dsMask, dsMesh, landMask, removable):
    print("Step 3: Perform flood fill, starting from open ocean.")

    # init flood fill to 0 for water, -1 for land, 1 for known open ocean
    floodFill = xarray.where(
        numpy.logical_and(removable, numpy.logical_not(landMask)), 0, -1)

    latCell = dsMesh.latCell
    lonCell = dsMesh.lonCell

    cellsOnCell = dsMesh.cellsOnCell

    # North Pole
    mask = latCell > 84.0
    openOceanMask = mask

    # Arctic
    mask = numpy.logical_and(
        numpy.logical_and(lonCell > 160.0, lonCell < 230.0),
        latCell > 73.0)
    openOceanMask = numpy.logical_or(openOceanMask, mask)

    # North Atlantic
    mask = numpy.logical_and(
        numpy.logical_and(lonCell > 315.0, lonCell < 340.0),
        numpy.logical_and(latCell > 15.0, latCell < 45.0))
    openOceanMask = numpy.logical_or(openOceanMask, mask)
    mask = numpy.logical_and(
        numpy.logical_and(lonCell > 290.0, lonCell < 300.0),
        numpy.logical_and(latCell > 72.0, latCell < 75.0))
    openOceanMask = numpy.logical_or(openOceanMask, mask)
    mask = numpy.logical_and(
        numpy.logical_and(lonCell > 0.0, lonCell < 10.0),
        numpy.logical_and(latCell > 70.0, latCell < 75.0))
    openOceanMask = numpy.logical_or(openOceanMask, mask)

    # North Pacific
    mask = numpy.logical_and(
        numpy.logical_and(lonCell > 150.0, lonCell < 225.0),
        numpy.logical_and(latCell > 0.0, latCell < 45.0))
    openOceanMask = numpy.logical_or(openOceanMask, mask)

    # South Atlantic
    mask = numpy.logical_and(
        numpy.logical_and(lonCell > 0.0, lonCell < 5.0),
        numpy.logical_and(latCell > -60.0, latCell < 0.0))
    openOceanMask = numpy.logical_or(openOceanMask, mask)

    # South Pacific
    mask = numpy.logical_and(
        numpy.logical_and(lonCell > 180.0, lonCell < 280.0),
        numpy.logical_and(latCell > -60.0, latCell < -10.0))
    openOceanMask = numpy.logical_or(openOceanMask, mask)

    # Southern Ocean
    mask = numpy.logical_and(
        numpy.logical_and(lonCell > 0.0, lonCell < 165.0),
        numpy.logical_and(latCell > -60.0, latCell < -45.0))
    openOceanMask = numpy.logical_or(openOceanMask, mask)

    mask = numpy.logical_and(floodFill == 0, openOceanMask)
    floodFill[mask] = 1

    nFloodableCells = numpy.count_nonzero(floodFill == 0)
    print("  Initial number of flood cells: {}".format(nFloodableCells))

    dsMask.landMaskDiagnostic[floodFill == 1] = 5

    # sweep over neighbors of known open ocean points
    for iSweep in range(dsMesh.sizes['nCells']):

        newFloodCellsThisSweep = 0
        mask = floodFill == 0
        cellIndices = numpy.nonzero(mask.values)[0]
        for iCellOnCell in range(cellsOnCell.shape[1]):
            neighbors = cellsOnCell[cellIndices, iCellOnCell]
            filledNeighbors = numpy.logical_and(neighbors >= 0,
                                                floodFill[neighbors] == 1)
            fillIndices = cellIndices[filledNeighbors.values]
            if(len(fillIndices) > 0):
                floodFill[fillIndices] = 1
                newFloodCellsThisSweep += len(fillIndices)

        print("  Sweep {} new flood cells this sweep: {}".format(
            iSweep, newFloodCellsThisSweep))

        if (newFloodCellsThisSweep == 0):
            break

    oceanMask = (floodFill == 1)

    print('oceanMask:', numpy.count_nonzero(oceanMask))

    return oceanMask


def _revert_cells_with_connected_edges(dsMask, dsMesh, oceanMask, landMask,
                                       removable, nSweeps):
    print("Step 4: Searching for land-locked cells, step 3: revert cells with "
          "connected active edges")

    cellsOnCell = dsMesh.cellsOnCell
    nextCellsOnCell = dsMesh.nextCellsOnCell
    prevCellsOnCell = dsMesh.prevCellsOnCell

    for iSweep in range(nSweeps):
        landMaskNew = numpy.array(landMask)

        # only remove a cell that was added in Step 2,
        # _remove_cells_with_isolated_edges2
        mask = numpy.logical_and(removable, dsMask.landMaskDiagnostic == 3)

        notLand = numpy.logical_not(landMask)
        valid = numpy.logical_and(mask, cellsOnCell >= 0)
        oceanEdge = numpy.logical_and(valid, oceanMask[cellsOnCell])
        valid = numpy.logical_and(mask, nextCellsOnCell >= 0)
        activeNextEdge = numpy.logical_and(valid, notLand[nextCellsOnCell])
        valid = numpy.logical_and(mask, prevCellsOnCell >= 0)
        activePrevEdge = numpy.logical_and(valid, notLand[prevCellsOnCell])

        reactivate = numpy.any(
            numpy.logical_and(
                oceanEdge,
                numpy.logical_or(activePrevEdge, activeNextEdge)), axis=1)

        landLockedCounter = numpy.count_nonzero(reactivate)
        if landLockedCounter > 0:
            landMaskNew[reactivate] = 0
            dsMask.regionCellMasks[reactivate, 0] = 0
            oceanMask[reactivate] = 1
            dsMask.landMaskDiagnostic[reactivate] = 4

        landMask = landMaskNew
        print("  Sweep: {} Number of land-locked cells returned: {}".format(
            iSweep + 1, landLockedCounter))
        if landLockedCounter == 0:
            break

    return landMask
