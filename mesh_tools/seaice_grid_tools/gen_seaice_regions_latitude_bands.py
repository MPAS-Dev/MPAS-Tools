#!/usr/bin/env python

import argparse
import numpy as np
from netCDF4 import Dataset
import math

radiansToDegrees = 180.0 / math.pi

# parsing
parser = argparse.ArgumentParser(description='Create partition regions based on latitude bands')

parser.add_argument('-m', '--mesh',      dest="meshFilename",   required=True,  help='MPAS mesh file')
parser.add_argument('-r', '--regions',   dest="regionFilename", required=True,  help='output region file')
parser.add_argument('-l', '--latitudes', dest="latitudes",      required=False, help='latitude band limits degrees not including north/south poles', nargs="+", type=float)

args = parser.parse_args()

# required arguments
meshFilename = args.meshFilename
regionFilename = args.regionFilename
if (args.latitudes is None):
    latitudes = np.array([-60.0,60.0])
else:
    latitudes = np.array(sorted(args.latitudes))


# load mesh file
mesh = Dataset(meshFilename,"r")
nCells = len(mesh.dimensions["nCells"])
latCell = mesh.variables["latCell"][:]
mesh.close()


# region variable
region = np.zeros(nCells,dtype="i")

nRegions = len(latitudes)+1

latitudeMin = np.zeros(nRegions)
latitudeMax = np.zeros(nRegions)

latitudeMin[0]  = -100.0
latitudeMax[0]  = latitudes[0]
latitudeMin[-1] = latitudes[-1]
latitudeMax[-1] =  100.0

for iRegion in range(1,nRegions-1):
    latitudeMin[iRegion] = latitudes[iRegion-1]
    latitudeMax[iRegion] = latitudes[iRegion]

print("nRegions:     ", nRegions)
print("latitudeMin:  ", latitudeMin)
print("latitudeMax:  ", latitudeMax)

for iRegion in range(0,nRegions):
    for iCell in range(0,nCells):
        if (latCell[iCell] * radiansToDegrees >= latitudeMin[iRegion] and
            latCell[iCell] * radiansToDegrees <  latitudeMax[iRegion]):
            region[iCell] = iRegion

# output regions
regionFile = Dataset(regionFilename,"w",format="NETCDF3_CLASSIC")

regionFile.nRegions = nRegions

regionFile.createDimension("nCells",nCells)

regionVar = regionFile.createVariable("region","i",dimensions=["nCells"])
regionVar[:] = region[:]

regionFile.close()
