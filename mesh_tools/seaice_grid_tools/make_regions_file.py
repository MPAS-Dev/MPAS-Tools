from __future__ import print_function
from netCDF4 import Dataset
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def make_regions_file(filenameIcePresent, filenameMesh, regionType, varname, limit, filenameOut):

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

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    # input parsing
    parser = argparse.ArgumentParser(description='Make a partition regions file.')

    parser.add_argument('-p', '--presence', dest="filenameIcePresent", required=True, help='ice presence file')
    parser.add_argument('-m', '--mesh',     dest="filenameMesh",       required=True, help='MPAS mesh file')
    # region type options: "two_region_eq", "three_region", "three_region_eq", "five_region_eq"
    parser.add_argument('-t', '--type',     dest="regionType",         required=True, help='region type')
    parser.add_argument('-v', '--varname',  dest="varname",            required=True, help='presence var name')
    parser.add_argument('-o', '--output',   dest="filenameOut",        required=True, help='output regions file')
    parser.add_argument('-l', '--limit',    dest="limit",              required=False, default=0.5, type=float, help='presence region limit')

    args = parser.parse_args()

    make_regions_file(args.filenameIcePresent, args.filenameMesh, args.regionType, args.varname, args.limit, args.filenameOut)
