#!/usr/bin/env python
"""
Name: add_land_locked_cells_to_mask.py
Author: Mark Petersen, Adrian Turner

Find ocean cells that are land-locked, and alter the cell
mask so that they are counted as land cells.
"""
import os
import shutil
from netCDF4 import Dataset
import numpy as np
import argparse

def removeFile(fileName):
    try:
        os.remove(fileName)
    except OSError:
        pass

parser = \
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--input_mask_file", dest="input_mask_filename",
                    help="Mask file that includes cell and edge masks.",
                    metavar="INPUTMASKFILE", required=True)
parser.add_argument("-o", "--output_mask_file", dest="output_mask_filename",
                    help="Mask file that includes cell and edge masks.",
                    metavar="OUTPUTMASKFILE", required=True)
parser.add_argument("-m", "--mesh_file", dest="mesh_filename",
                    help="MPAS Mesh filename.", metavar="MESHFILE",
                    required=True)
parser.add_argument("-l", "--latitude_threshold", dest="latitude_threshold",
                    help="Minimum latitude, in degrees, for transect widening.",
                    required=False, type=float, default=43.0)
parser.add_argument("-n", "--number_sweeps", dest="nSweeps",
                    help="Maximum number of sweeps to search for land-locked cells.",
                    required=False, type=int, default=10)
args = parser.parse_args()

latitude_threshold_radians = args.latitude_threshold*3.1415/180.

# Obtain mesh variables
meshFile = Dataset(args.mesh_filename, "r")
nCells = len(meshFile.dimensions["nCells"])
maxEdges = len(meshFile.dimensions["maxEdges"])
cellsOnCell = meshFile.variables["cellsOnCell"][:, :]
nEdgesOnCell = meshFile.variables["nEdgesOnCell"][:]
latCell = meshFile.variables["latCell"][:]
lonCell = meshFile.variables["lonCell"][:]
meshFile.close()

removeFile(args.output_mask_filename)
shutil.copyfile(args.input_mask_filename,args.output_mask_filename)

# Obtain original cell mask from input file
inputMaskFile = Dataset(args.input_mask_filename, "r")
nRegions = len(inputMaskFile.dimensions["nRegions"])
regionCellMasks = inputMaskFile.variables["regionCellMasks"][:, :]
# set landMask = flattened regionCellMasks
landMask = np.amax(regionCellMasks, axis=1)
inputMaskFile.close()

# Open output file
outputMaskFile = Dataset(args.output_mask_filename, "a")
landMaskDiagnostic = outputMaskFile.createVariable("landMaskDiagnostic", "i", dimensions=("nCells"))

print "Running add_land_locked_cells_to_mask.py.  Total number of cells: ", nCells

# use np.array, as simple = makes a pointer
landMaskNew = np.array(landMask)
activeEdgeSum = np.zeros(maxEdges, dtype="i")

# Removable cells are ocean cells outside of latitude threshold
removableCellIndex = np.zeros(nCells, dtype="i")
nRemovableCells = 0

print "Step 1: Searching for land-locked cells.  Remove cells that only have isolated active edges."
landLockedCounter = 0
for iCell in range(nCells):
    landMaskDiagnostic[iCell] = landMask[iCell]
    # skip if outside latitude threshold or if this is already a land cell
    if abs(latCell[iCell]) < latitude_threshold_radians or landMask[iCell] == 1:
        continue
    removableCellIndex[nRemovableCells] = iCell
    nRemovableCells += 1
    activeEdgeSum[:] = 0
    for iEdgeOnCell in range(nEdgesOnCell[iCell]):
        # check if neighbor is an ocean cell (landMask=0)
        # subtract 1 to convert 1-base to 0-base:
        if landMask[cellsOnCell[iCell, iEdgeOnCell]-1] == 0:
            activeEdgeSum[iEdgeOnCell] += 1
            # % is modulo operator:
            iP1 = (iEdgeOnCell + 1) % nEdgesOnCell[iCell]
            activeEdgeSum[iP1] += 1

    if np.amax(activeEdgeSum[0:nEdgesOnCell[iCell]]) == 1:
        outputMaskFile['regionCellMasks'][iCell, 1] = 1
        landLockedCounter += 1
        landMaskNew[iCell] = 1
        landMaskDiagnostic[iCell] = 2

landMask[:] = landMaskNew[:]
print "  Number of landLocked cells: ", landLockedCounter

print "Step 2: Searching for land-locked cells. Remove cells that have any isolated active edges."
for iSweep in range(args.nSweeps):
    landLockedCounter = 0
    for iRemovableCell in range(0, nRemovableCells):
        iCell = removableCellIndex[iRemovableCell]
        if landMask[iCell] == 1:
            continue
        for iEdgeOnCell in range(nEdgesOnCell[iCell]):
            # check if neighbor is an ocean cell (landMask=0)
            # subtract 1 to convert 1-base to 0-base:
            if landMask[cellsOnCell[iCell, iEdgeOnCell]-1] == 0:
                # % is modulo operator:
                iP1 = (iEdgeOnCell + 1) % nEdgesOnCell[iCell]
                iM1 = (iEdgeOnCell - 1) % nEdgesOnCell[iCell]
                # Is this neighbor's two neighbors to left and right land?
                # if so, sum of masks is two.
                # subtract 1 to convert 1-base to 0-base:
                if (landMask[cellsOnCell[iCell, iP1]-1]
                       + landMask[cellsOnCell[iCell, iM1]-1]) == 2:
                    landLockedCounter += 1
                    landMaskNew[iCell] = 1
                    outputMaskFile['regionCellMasks'][iCell, 1] = 1
                    landMaskDiagnostic[iCell] = 3
                    # once we remove this cell, we can quit checking over edges
                    break

    landMask[:] = landMaskNew[:]
    print "  Sweep: ", iSweep+1, "Number of landLocked cells removed: ", landLockedCounter
    if landLockedCounter == 0:
        break

print "Step 3: Perform flood fill, starting from open ocean."
floodFill = np.zeros(nCells, dtype="i")
floodableCellIndex = np.zeros(nCells, dtype="i")
nFloodableCells = 0
floodFill[:] = -1
d2r = 3.1415/180.0

# init flood fill to 0 for water, -1 for land, 1 for known open ocean regions
for iRemovableCell in range(0, nRemovableCells):
    iCell = removableCellIndex[iRemovableCell]
    if (landMaskDiagnostic[iCell] == 0):
        floodFill[iCell] = 0
        if (latCell[iCell] > 84.0*d2r  # North Pole
            or lonCell[iCell] > 160.0*d2r and lonCell[iCell] < 230.0*d2r and latCell[iCell] >  73.0*d2r   # Arctic
            or lonCell[iCell] > 315.0*d2r and lonCell[iCell] < 340.0*d2r and latCell[iCell] >  15.0*d2r and latCell[iCell] <  45.0*d2r  # North Atlantic
            or lonCell[iCell] > 290.0*d2r and lonCell[iCell] < 300.0*d2r and latCell[iCell] >  72.0*d2r and latCell[iCell] <  75.0*d2r  # North Atlantic
            or lonCell[iCell] >   0.0*d2r and lonCell[iCell] <  10.0*d2r and latCell[iCell] >  70.0*d2r and latCell[iCell] <  75.0*d2r  # North Atlantic 2
            or lonCell[iCell] > 150.0*d2r and lonCell[iCell] < 225.0*d2r and latCell[iCell] >   0.0*d2r and latCell[iCell] <  45.0*d2r  # North Pacific
            or lonCell[iCell] >   0.0*d2r and lonCell[iCell] <   5.0*d2r and latCell[iCell] > -60.0*d2r and latCell[iCell] <   0.0*d2r  # South Atlantic
            or lonCell[iCell] > 180.0*d2r and lonCell[iCell] < 280.0*d2r and latCell[iCell] > -60.0*d2r and latCell[iCell] < -10.0*d2r  # South Pacific
            or lonCell[iCell] >   0.0*d2r and lonCell[iCell] < 165.0*d2r and latCell[iCell] > -60.0*d2r and latCell[iCell] < -45.0*d2r):  # Southern Ocean
            floodFill[iCell] = 1
            landMaskDiagnostic[iCell] = 5  # indicates seed region
        else:
            floodableCellIndex[nFloodableCells] = iCell
            nFloodableCells += 1
print "  Initial number of flood cells: ", nFloodableCells

# sweep over neighbors of known open ocean points
for iSweep in range(0, nCells):
    newFloodCellsThisSweep = 0

    for iFloodableCell in range(0, nFloodableCells):
        iCell = floodableCellIndex[iFloodableCell]
        if (floodFill[iCell] == 0):

            for iCellOnCellSweep in range(0, nEdgesOnCell[iCell]):
                iCellNeighbor = cellsOnCell[iCell, iCellOnCellSweep]-1

                if (floodFill[iCellNeighbor] == 1):
                    floodFill[iCell] = 1
                    newFloodCellsThisSweep += 1
                    break

    print "  Sweep ", iSweep, " new flood cells this sweep: ", newFloodCellsThisSweep

    if (newFloodCellsThisSweep == 0):
        break

oceanMask = np.zeros(nCells, dtype="i")
for iCell in range(0, nCells):
    if (floodFill[iCell] == 1):
        oceanMask[iCell] = 1

print "Step 4: Searching for land-locked cells, step 3: revert cells with connected active edges"
for iSweep in range(args.nSweeps):
    landLockedCounter = 0
    for iRemovableCell in range(0, nRemovableCells):
        iCell = removableCellIndex[iRemovableCell]
        # only remove a cell that was added in lats round (red cells)
        if landMaskDiagnostic[iCell] == 3:
            for iEdgeOnCell in range(nEdgesOnCell[iCell]):
                # check if neighbor is an ocean cell (landMask=0)
                # subtract 1 to convert 1-base to 0-base:
                if oceanMask[cellsOnCell[iCell, iEdgeOnCell]-1] == 1:
                    # % is modulo operator:
                    iP1 = (iEdgeOnCell + 1) % nEdgesOnCell[iCell]
                    iM1 = (iEdgeOnCell - 1) % nEdgesOnCell[iCell]
                    # Is either of this neighbor's two neighbors to left and right ocean?
                    # if so, sum of masks is two.
                    # subtract 1 to convert 1-base to 0-base:
                    if (landMask[cellsOnCell[iCell, iP1]-1] == 0
                            or landMask[cellsOnCell[iCell, iM1]-1] == 0):
                        landLockedCounter += 1
                        landMaskNew[iCell] = 0
                        outputMaskFile['regionCellMasks'][iCell, 1] = 0
                        landMaskDiagnostic[iCell] = 4
                        oceanMask[iCell] = 1
                        # once we remove this cell, we can quit checking over edges
                        break

    landMask[:] = landMaskNew[:]
    print "  Sweep: ", iSweep+1, "Number of land-locked cells returned: ", landLockedCounter
    if landLockedCounter == 0:
        break

outputMaskFile.close()
