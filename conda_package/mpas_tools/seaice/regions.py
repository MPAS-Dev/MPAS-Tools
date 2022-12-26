from __future__ import print_function
from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def make_regions_file(filenameIcePresent, filenameMesh, regionType, varname, limit, filenameOut):
    """"

    Parameters
    ----------
    filenameIcePresent : str
        A filename for a file containing the field specified in ``varname``
        that determines where ice may be present

    filenameMesh : str
        The name of a file containing the MPAS-Seaice mesh

    regionType : {"three_region", "two_region_eq", "three_region_eq", "five_region_eq"}
        The type of regions to write

    varname : str
        The name of the variable that determines where ice might be present

    limit : float
        For ``regionType`` either ``three_region`` or ``five_region_eq``,
        the  value of the ``varname`` field above which ice is considered to
        be present.  For other ``regionType`` values, the limit is always 0.5
        and this parameter is ignored.

    filenameOut : str
        The NetCDF filename to write the resulting ``region`` field to
    """

    # load ice presence
    print(filenameIcePresent)
    fileIcePresent = Dataset(filenameIcePresent,"r")

    icePresence = fileIcePresent.variables[varname][:]

    fileIcePresent.close()

    # load mesh data
    fileMesh = Dataset(filenameMesh,"r")

    nCells = len(fileMesh.dimensions["nCells"])

    latCell = fileMesh.variables["latCell"][:]

    fileMesh.close()

    # output region mask
    region = np.zeros(nCells,dtype="i")


    # no equatorial removal
    if (regionType == "three_region"):

        nRegions = 3

        for iCell in range(0,nCells):

            if   (icePresence[iCell] > limit and latCell[iCell] < 0.0):
                region[iCell] = 0
            elif (icePresence[iCell] > limit and latCell[iCell] >= 0.0):
                region[iCell] = 2
            else:
                region[iCell] = 1

    elif (regionType == "two_region_eq"):

        nRegions = 2

        for iCell in range(0,nCells):

            if   (icePresence[iCell] > 0.5 and latCell[iCell] < 0.0):
                region[iCell] = 0
            elif (icePresence[iCell] > 0.5 and latCell[iCell] >= 0.0):
                region[iCell] = 0
            else:
                region[iCell] = 1

    elif (regionType == "three_region_eq"):

        nRegions = 3

        for iCell in range(0,nCells):

            if   (icePresence[iCell] > 0.5 and latCell[iCell] < 0.0):
                region[iCell] = 0
            elif (icePresence[iCell] > 0.5 and latCell[iCell] >= 0.0):
                region[iCell] = 2
            else:
                region[iCell] = 1

    elif (regionType == "five_region_eq"):

        nRegions = 5

        for iCell in range(0,nCells):

            if   (icePresence[iCell] > limit and latCell[iCell] >= 0.0):
                region[iCell] = 4
            elif (icePresence[iCell] < limit and icePresence[iCell] > 0.5 and latCell[iCell] >= 0.0):
                region[iCell] = 3
            elif (icePresence[iCell] > limit and latCell[iCell] < 0.0):
                region[iCell] = 0
            elif (icePresence[iCell] < limit and icePresence[iCell] > 0.5 and latCell[iCell] < 0.0):
                region[iCell] = 1
            else:
                region[iCell] = 2


    # output
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.nRegions = nRegions

    fileOut.createDimension("nCells",nCells)

    regionVar = fileOut.createVariable("region","i",dimensions=["nCells"])
    regionVar[:] = region[:]

    fileOut.close()
