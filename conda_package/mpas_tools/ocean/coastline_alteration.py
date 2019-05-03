from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
from netCDF4 import Dataset
import os
import shutil


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


def add_land_locked_cells_to_mask(input_mask_filename, output_mask_filename,
                                  mesh_filename, latitude_threshold=43.0,
                                  nSweeps=10):
    '''
    Find ocean cells that are land-locked, and alter the cell mask so that they
    are counted as land cells.

    Parameters
    ----------
    input_mask_filename : str
        Mask file that includes cell and edge masks.

    output_mask_filename : str
        Mask file that includes cell and edge masks.

    mesh_filename : str
        MPAS Mesh filename.

    latitude_threshold : float, optional
        Minimum latitude, in degrees, for transect widening.

    nSweeps : int, optional
        Maximum number of sweeps to search for land-locked cells.
    '''

    # Obtain mesh variables
    meshFile = Dataset(mesh_filename, "r")
    nCells = len(meshFile.dimensions["nCells"])
    cellsOnCell = meshFile.variables["cellsOnCell"][:, :] - 1
    nEdgesOnCell = meshFile.variables["nEdgesOnCell"][:]
    latCell = numpy.rad2deg(meshFile.variables["latCell"][:])
    lonCell = numpy.rad2deg(meshFile.variables["lonCell"][:])
    meshFile.close()

    _remove_file(output_mask_filename)
    shutil.copyfile(input_mask_filename, output_mask_filename)

    # Obtain original cell mask from input file
    inputMaskFile = Dataset(input_mask_filename, "r")
    regionCellMasks = inputMaskFile.variables["regionCellMasks"][:, :]
    # set landMask = flattened regionCellMasks
    landMask = numpy.amax(regionCellMasks, axis=1)
    inputMaskFile.close()

    # Open output file
    outputMaskFile = Dataset(output_mask_filename, "a")
    landMaskDiagnostic = outputMaskFile.createVariable(
        "landMaskDiagnostic", "i", dimensions=("nCells"))

    regionCellMasks = outputMaskFile['regionCellMasks']

    print("Running add_land_locked_cells_to_mask.py.  Total number of cells: "
          "{}".format(nCells))

    landMask, removable = _remove_cells_with_isolated_edges1(
        landMask, landMaskDiagnostic, regionCellMasks, latCell, nEdgesOnCell,
        cellsOnCell, nCells, latitude_threshold)
    landMask = _remove_cells_with_isolated_edges2(
        landMask, landMaskDiagnostic, regionCellMasks, removable,
        nEdgesOnCell, cellsOnCell, nCells, nSweeps)
    oceanMask = _flood_fill(landMask, landMaskDiagnostic, removable, lonCell,
                            latCell, nEdgesOnCell, cellsOnCell, nCells)
    landMask = _revert_cells_with_connected_edges(
        oceanMask, landMask, landMaskDiagnostic, regionCellMasks, removable,
        nEdgesOnCell, cellsOnCell, nCells, nSweeps)
    outputMaskFile.close()


def _remove_cells_with_isolated_edges1(landMask, landMaskDiagnostic,
                                       regionCellMasks, latCell, nEdgesOnCell,
                                       cellsOnCell, nCells,
                                       latitude_threshold):
    print("Step 1: Searching for land-locked cells.  Remove cells that only "
          "have isolated active edges.")

    landMaskNew = numpy.array(landMask)

    landMaskDiagnostic[:] = landMask

    removable = numpy.logical_and(numpy.abs(latCell) >= latitude_threshold,
                                  landMask == 0)

    nextCellsOnCell = numpy.zeros(cellsOnCell.shape, int)

    for iEdgeOnCell in range(nextCellsOnCell.shape[1]):
        iP1 = numpy.mod(iEdgeOnCell + 1, nEdgesOnCell)
        nextCellsOnCell[:, iEdgeOnCell] = \
            cellsOnCell[numpy.arange(nCells), iP1]

    valid = numpy.logical_and(removable.reshape(nCells, 1),
                              cellsOnCell >= 0)

    active = numpy.logical_not(landMask)
    activeEdge = numpy.logical_and(valid, active[cellsOnCell])
    activeNextEdge = numpy.logical_and(valid, active[nextCellsOnCell])

    # which vertices have adjacent active edges on this cell?
    activeAdjacentEdges = numpy.logical_and(activeEdge, activeNextEdge)

    # which removable cells have no pairs of adjacent active cells?
    noActiveAdjacentEdges = numpy.logical_and(
        removable, numpy.logical_not(numpy.any(activeAdjacentEdges, axis=1)))

    landMaskNew[noActiveAdjacentEdges] = 1
    landLockedCounter = numpy.count_nonzero(noActiveAdjacentEdges)

    regionCellMasks[:, 0] = numpy.maximum(regionCellMasks[:, 0],
                                          noActiveAdjacentEdges)

    landMaskDiagnostic[noActiveAdjacentEdges] = 2

    print("  Number of landLocked cells: {}".format(landLockedCounter))

    return landMaskNew, removable


def _remove_cells_with_isolated_edges2(landMask, landMaskDiagnostic,
                                       regionCellMasks, removable,
                                       nEdgesOnCell, cellsOnCell, nCells,
                                       nSweeps):
    print("Step 2: Searching for land-locked cells. Remove cells that have "
          "any isolated active edges.")

    nextCellsOnCell = numpy.zeros(cellsOnCell.shape, int)
    prevCellsOnCell = numpy.zeros(cellsOnCell.shape, int)

    for iEdgeOnCell in range(nextCellsOnCell.shape[1]):
        iP1 = numpy.mod(iEdgeOnCell + 1, nEdgesOnCell)
        nextCellsOnCell[:, iEdgeOnCell] = \
            cellsOnCell[numpy.arange(nCells), iP1]
        iM1 = numpy.mod(iEdgeOnCell - 1, nEdgesOnCell)
        prevCellsOnCell[:, iEdgeOnCell] = \
            cellsOnCell[numpy.arange(nCells), iM1]

    for iSweep in range(nSweeps):
        landLockedCounter = 0
        landMaskNew = numpy.array(landMask)

        mask = numpy.logical_and(removable,
                                 landMask == 0)

        active = numpy.logical_not(landMask)
        valid = numpy.logical_and(mask.reshape(nCells, 1),
                                  cellsOnCell >= 0)
        activeEdge = numpy.logical_and(valid, active[cellsOnCell])
        valid = numpy.logical_and(mask.reshape(nCells, 1),
                                  nextCellsOnCell >= 0)
        activeNextEdge = numpy.logical_and(valid, active[nextCellsOnCell])
        valid = numpy.logical_and(mask.reshape(nCells, 1),
                                  prevCellsOnCell >= 0)
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
            regionCellMasks[landLockedCells, 0] = 1
            landMaskDiagnostic[landLockedCells] = 3

        landMask = landMaskNew
        print("  Sweep: {} Number of landLocked cells removed: {}".format(
            iSweep + 1, landLockedCounter))
        if landLockedCounter == 0:
            break

    return landMask


def _flood_fill(landMask, landMaskDiagnostic, removable, lonCell, latCell,
                nEdgesOnCell, cellsOnCell, nCells):
    print("Step 3: Perform flood fill, starting from open ocean.")

    # init flood fill to 0 for water, -1 for land, 1 for known open ocean
    floodFill = -1*numpy.ones(nCells, dtype="i")
    mask = numpy.logical_and(removable, landMask == 0)
    floodFill[mask] = 0

    openOceanMask = numpy.zeros(nCells, bool)

    # North Pole
    mask = latCell > 84.0
    openOceanMask = numpy.logical_or(openOceanMask, mask)

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

    landMaskDiagnostic[floodFill == 1] = 5

    # sweep over neighbors of known open ocean points
    for iSweep in range(0, nCells):

        newFloodCellsThisSweep = 0
        mask = floodFill == 0
        for iCellOnCell in range(cellsOnCell.shape[1]):
            neighbors = cellsOnCell[:, iCellOnCell]
            fill = numpy.logical_and(
                mask,
                numpy.logical_and(neighbors >= 0, floodFill[neighbors] == 1))
            floodFill[fill] = 1
            newFloodCellsThisSweep += numpy.count_nonzero(fill)

        print("  Sweep {} new flood cells this sweep: {}".format(
            iSweep, newFloodCellsThisSweep))

        if (newFloodCellsThisSweep == 0):
            break

    oceanMask = (floodFill == 1)

    print('oceanMask:', numpy.count_nonzero(oceanMask))

    return oceanMask


def _revert_cells_with_connected_edges(oceanMask, landMask, landMaskDiagnostic,
                                       regionCellMasks, removable,
                                       nEdgesOnCell, cellsOnCell, nCells,
                                       nSweeps):
    print("Step 4: Searching for land-locked cells, step 3: revert cells with "
          "connected active edges")

    nextCellsOnCell = numpy.zeros(cellsOnCell.shape, int)
    prevCellsOnCell = numpy.zeros(cellsOnCell.shape, int)

    for iEdgeOnCell in range(nextCellsOnCell.shape[1]):
        iP1 = numpy.mod(iEdgeOnCell + 1, nEdgesOnCell)
        nextCellsOnCell[:, iEdgeOnCell] = \
            cellsOnCell[numpy.arange(nCells), iP1]
        iM1 = numpy.mod(iEdgeOnCell - 1, nEdgesOnCell)
        prevCellsOnCell[:, iEdgeOnCell] = \
            cellsOnCell[numpy.arange(nCells), iM1]

    for iSweep in range(nSweeps):
        landMaskNew = numpy.array(landMask)

        # only remove a cell that was added in Step 2,
        # _remove_cells_with_isolated_edges2
        mask = numpy.logical_and(removable, landMaskDiagnostic[:] == 3)

        valid = numpy.logical_and(mask.reshape(nCells, 1),
                                  cellsOnCell >= 0)
        oceanEdge = numpy.logical_and(valid, oceanMask[cellsOnCell])
        valid = numpy.logical_and(mask.reshape(nCells, 1),
                                  nextCellsOnCell >= 0)
        activeNextEdge = numpy.logical_and(valid,
                                           landMask[nextCellsOnCell] == 0)
        valid = numpy.logical_and(mask.reshape(nCells, 1),
                                  prevCellsOnCell >= 0)
        activePrevEdge = numpy.logical_and(valid,
                                           landMask[prevCellsOnCell] == 0)

        reactivate = numpy.any(
            numpy.logical_and(
                oceanEdge,
                numpy.logical_or(activePrevEdge, activeNextEdge)), axis=1)

        landLockedCounter = numpy.count_nonzero(reactivate)
        if landLockedCounter > 0:
            landMaskNew[reactivate] = 0
            regionCellMasks[reactivate, 0] = 0
            oceanMask[reactivate] = 1
            landMaskDiagnostic[reactivate] = 4

        landMask = landMaskNew
        print("  Sweep: {} Number of land-locked cells returned: {}".format(
            iSweep + 1, landLockedCounter))
        if landLockedCounter == 0:
            break

    return landMask


def _remove_file(fileName):
    try:
        os.remove(fileName)
    except OSError:
        pass
