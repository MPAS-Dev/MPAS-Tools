#!/usr/bin/env python
# Create a SCRIP file from a planar rectanfular mesh.
# See for details: http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_5_2_0rp1/ESMF_refdoc/node3.html#SECTION03024000000000000000

import sys
import netCDF4
import numpy as np
from optparse import OptionParser


print ("== Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser()
parser.description = "This script takes an MPAS grid file and generates a SCRIP grid file."
parser.add_option("-i", "--input", dest="inputFile", help="input grid file name used as input.", default="input.nc", metavar="FILENAME")
parser.add_option("-s", "--scrip", dest="scripFile", help="SCRIP grid file to output.", default="scrip.nc", metavar="FILENAME")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

if not options.inputFile:
    sys.exit('Error: Data input grid file is required.  Specify with -c command line argument.')
if not options.scripFile:
    sys.exit('Error: SCRIP output grid file is required.  Specify with -s command line argument.')
print ('') # make a space in stdout before further output

# ===================================

fin = netCDF4.Dataset(options.inputFile, 'r')
fout = netCDF4.Dataset(options.scripFile, 'w')  # This will clobber existing files

# Get info from input file
nx = len(fin.dimensions['x'])
ny = len(fin.dimensions['y'])

# Write to output file
# Dimensions
fout.createDimension("grid_size", nx * ny)
fout.createDimension("grid_corners", 4 )

print('grid rank is 2')
fout.createDimension("grid_rank", 2)

# Variables
grid_center_lat = fout.createVariable('grid_center_lat', 'f8', ('grid_size',))
grid_center_lat.units = 'degrees'
grid_center_lon = fout.createVariable('grid_center_lon', 'f8', ('grid_size',))
grid_center_lon.units = 'degrees'
grid_corner_lat = fout.createVariable('grid_corner_lat', 'f8', ('grid_size', 'grid_corners'))
grid_corner_lat.units = 'degrees'
grid_corner_lon = fout.createVariable('grid_corner_lon', 'f8', ('grid_size', 'grid_corners'))
grid_corner_lon.units = 'degrees'
grid_imask = fout.createVariable('grid_imask', 'i4', ('grid_size',))
grid_imask.units = 'unitless'
grid_dims = fout.createVariable('grid_dims', 'i4', ('grid_rank',))


grid_center_lat[:] = fin.variables['lat'][:].flatten()
grid_center_lon[:] = fin.variables['lon'][:].flatten()

lat_bnds = fin.variables['lat_bnds'][:]
grid_corner_lat[:] = lat_bnds.reshape(-1, lat_bnds.shape[-1])
lon_bnds = fin.variables['lon_bnds'][:]
grid_corner_lon[:] = lon_bnds.reshape(-1, lon_bnds.shape[-1])

grid_imask[:] = 1  # For now, assume we don't want to mask anything out - but eventually may want to exclude certain cells from the input mesh during interpolation

# set the grid dimension based on the grid rank
grid_dims[:] = [nx , ny]

fin.close()
fout.close()
print('scrip file generation complete')
