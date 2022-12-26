from netCDF4 import Dataset
import numpy as np
import math

from mpas_tools.cime.constants import constants


#-------------------------------------------------------------------------------

def extend_seaice_mask(filenameMesh,filenamePresence,extendDistance,unitSphere=False):
    """
    Add a field ``icePresenceExtended`` to ``filenamePresence`` if it doesn't
    already exist.  This field is the ``icePresence`` field extended by
    a distance of ``extendDistance``.

    Parameters
    ----------
    filenameMesh : str
        The filename of the MPAS-Seaice mesh

    filenamePresence : str
        A filename for a file containing an ``icePresence`` field to be
        extended and to which a new ``icePresenceExtended`` will be added

    extendDistance : float
        The distance in km to expand (no expansion is performed if
        ``extendDistance=0.0``)

    unitSphere : bool, optional
        Whether the mesh is provided on the unit sphere, rather than one with
        the radius of the Earth
    """

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

    if "icePresenceExtended" in filePresence.variables:
        # we're already done, probably from a previous call
        filePresence.close()
        return

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

    extendDistance = extendDistance * 1000.0 # convert from km to m

    if (unitSphere):
        earthRadius = constants['SHR_CONST_REARTH']
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
