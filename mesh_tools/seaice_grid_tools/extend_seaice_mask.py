from __future__ import print_function
from netCDF4 import Dataset
import numpy as np
import math
import argparse

#-------------------------------------------------------------------------------

def extend_seaice_mask(filenameMesh,filenamePresence,extendDistance,unitSphere=False):

    # mesh
    print("Load mesh...")
    fileMesh = Dataset(filenameMesh,"r")

    nCells = len(fileMesh.dimensions["nCells"])

    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    cellsOnCell = fileMesh.variables["cellsOnCell"][:]

    cellsOnCell[:] = cellsOnCell[:] - 1

    xCell = fileMesh.variables["xCell"][:]
    yCell = fileMesh.variables["yCell"][:]
    zCell = fileMesh.variables["zCell"][:]

    fileMesh.close()

    # presence
    print("Load ice presence...")
    filePresence = Dataset(filenamePresence,"r")

    icePresence = filePresence.variables["icePresence"][:]

    filePresence.close()

    # ice edge cells
    print("Get ice edge cells...")
    iceEdgeCell = np.zeros(nCells,dtype="i")

    for iCell in range(0,nCells):

        #if (iCell % 100000 == 0):
        #    print(iCell, " of ", nCells, " cells...")

        for iCellOnCell in range(0,nEdgesOnCell[iCell]):

            iCell2 = cellsOnCell[iCell,iCellOnCell]

            if (icePresence[iCell] == 1 and icePresence[iCell2] == 0):
                iceEdgeCell[iCell] = 1

    # only edge cells
    nEdgeCells = np.sum(iceEdgeCell)
    print("nEdgeCells: ", nEdgeCells)

    print("Get edge cell vector...")
    iCellEdge = np.zeros(nEdgeCells,dtype="i")

    iEdgeCell = 0
    for iCell in range(0, nCells):
        if (iceEdgeCell[iCell] == 1):
            iCellEdge[iEdgeCell] = iCell
            iEdgeCell = iEdgeCell + 1

    del iceEdgeCell


    # get edge positions
    print("Get edge positions...")
    xCellEdgeCell = np.zeros(nEdgeCells)
    yCellEdgeCell = np.zeros(nEdgeCells)
    zCellEdgeCell = np.zeros(nEdgeCells)

    for iEdgeCell in range(0,nEdgeCells):
        iCell = iCellEdge[iEdgeCell]
        xCellEdgeCell[iEdgeCell] = xCell[iCell]
        yCellEdgeCell[iEdgeCell] = yCell[iCell]
        zCellEdgeCell[iEdgeCell] = zCell[iCell]

    # find extended mask
    print("Find new extended mask...")

    earthRadius = 6371229.0
    extendDistance = extendDistance * 1000.0 # convert from km to m

    if (unitSphere):
        extendDistance = extendDistance / earthRadius

    distanceLimit = math.pow(extendDistance,2)

    icePresenceNew = icePresence.copy()

    if (extendDistance > 0.0):

        distances = np.zeros(nEdgeCells)

        for iCell in range(0,nCells):

            if (iCell % 100000 == 0):
                print(iCell, " of ", nCells, " cells...")

            distances = np.multiply((xCell[iCell] - xCellEdgeCell),(xCell[iCell] - xCellEdgeCell)) + \
                        np.multiply((yCell[iCell] - yCellEdgeCell),(yCell[iCell] - yCellEdgeCell)) + \
                        np.multiply((zCell[iCell] - zCellEdgeCell),(zCell[iCell] - zCellEdgeCell))

            if (distances[distances.argmin()] <= distanceLimit):
                icePresenceNew[iCell] = 1

    print("Load ice presence...")
    filePresence = Dataset(filenamePresence,"a")

    icePresenceVar = filePresence.createVariable("icePresenceExtended","d",dimensions=["nCells"])
    icePresenceVar[:] = icePresenceNew[:]

    filePresence.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extend the ice presence variable')

    parser.add_argument('-m', '--inputmesh',  dest="filenameMesh",     required=True,  help='MPAS mesh file for source regridding mesh')
    parser.add_argument('-p', '--presence',   dest="filenamePresence", required=True,  help='File with ice presence')
    parser.add_argument('-e', '--extenddist', dest="extendDistance",   required=True,  help='distance (km) to extend ice present region', type=float)
    parser.add_argument('-u', '--unitsphere', dest="unitSphere",       required=False, help='Is the mesh file a unit sphere', action='store_false')

    args = parser.parse_args()

    extend_seaice_mask(args.filenameMesh,args.filenamePresence,args.extendDistance,args.unitSphere)
