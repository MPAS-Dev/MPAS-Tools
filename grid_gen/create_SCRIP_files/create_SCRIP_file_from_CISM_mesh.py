#!/usr/bin/env python
# Create a SCRIP file from a CISM mesh.
# See for details: http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_5_2_0rp1/ESMF_refdoc/node3.html#SECTION03024000000000000000

import sys
import netCDF4
import numpy as np
from optparse import OptionParser
import pyproj

# ======== DEFINE PROJECTIONS =============
# Create empty dictionary to store projection definitions:
projections = dict()
# add more as needed:

# CISM's projection is as follows, with the vertical datum as EIGEN-GL04C geoid. 
# datum is actually EIGEN-GL04C but that is not an option in Proj.  Therefore using EGM08 which should be within ~1m everywhere (and 10-20 cm in most places)
# NOTE!!!!!!  egm08_25.gtx can be downloaded from:  http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx  and the path in the projection specification line should point to it!
#projections['gis-bamber'] = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids=./egm08_25.gtx')
projections['gis-bamber'] = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +ellps=WGS84')  # This version ignores the vertical datum shift, which should be a very small error for horizontal-only positions

# GIMP projection: This is also polar stereographic but with different standard parallel and using the WGS84 ellipsoid.
projections['gis-gimp'] = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

# BEDMAP2 projection
projections['ais-bedmap2'] = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')  # Note: BEDMAP2 elevations use EIGEN-GL04C geoid

# Standard Lat/Long
projections['latlon'] = pyproj.Proj(proj='latlong', datum='WGS84')
# ===================================



print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.description = "This script takes an MPAS grid file and generates a SCRIP grid file."
parser.add_option("-c", "--cism", dest="cismFile", help="CISM grid file name used as input.", default="cism.nc", metavar="FILENAME")
parser.add_option("-s", "--scrip", dest="scripFile", help="SCRIP grid file to output.", default="scrip.nc", metavar="FILENAME")
parser.add_option("-p", "--proj", dest="projection", help="projection used by the CISM file. Valid options are:"  + str(projections.keys()), metavar="PROJ")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

if not options.cismFile:
    sys.exit('Error: CISM input grid file is required.  Specify with -c command line argument.')
if not options.scripFile:
    sys.exit('Error: SCRIP output grid file is required.  Specify with -s command line argument.')
if not options.projection:
    sys.exit('Error: data projection required with -p or --proj command line argument. Valid options are: ' + str(projections.keys()))

print '' # make a space in stdout before further output


# ===============================================

fin = netCDF4.Dataset(options.cismFile, 'r')
fout = netCDF4.Dataset(options.scripFile, 'w')  # This will clobber existing files

# Get info from input file
x1 = fin.variables['x1'][:]
y1 = fin.variables['y1'][:]
nx = x1.size
ny = y1.size
dx = x1[1] - x1[0]
dy = y1[1] - y1[0]

# Write to output file
# Dimensions
fout.createDimension("grid_size", nx * ny)
fout.createDimension("grid_corners", 4 )
fout.createDimension("grid_rank", 1)

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


# Create matrices of x,y
print 'Building matrix version of x, y locations.'
x1matrix, y1matrix = np.meshgrid(x1, y1)
xc = np.append(x1[:] - dx/2.0, x1[-1] + dx/2.0 )  # get a copy of x1 that is on the staggered grid and includes both bounding edges
yc = np.append(y1[:] - dy/2.0, y1[-1] + dy/2.0 )  # get a copy of x1 that is on the staggered grid and includes both bounding edges
xcmatrix, ycmatrix = np.meshgrid(xc, yc)

# Unproject to lat/long for grid centers and grid corners
print 'Unprojecting.'

x1matrix_flat = x1matrix.flatten(order='C')  # Flatten using C indexing
y1matrix_flat = y1matrix.flatten(order='C')
grid_center_lon[:], grid_center_lat[:] = pyproj.transform(projections[options.projection], projections['latlon'], x1matrix_flat, y1matrix_flat, radians=False)


# Now fill in the corners in the right locations
stag_lon, stag_lat = pyproj.transform(projections[options.projection], projections['latlon'], xcmatrix, ycmatrix, radians=False)
print 'Filling in corners of each cell.'
grid_corner_lon_local = np.zeros( (nx * ny, 4) )  # It is WAYYY faster to fill in the array entry-by-entry in memory than to disk.
grid_corner_lat_local = np.zeros( (nx * ny, 4) )
for j in range(ny):
  for i in range(nx):
    iCell = j*nx + i

    grid_corner_lon_local[iCell, 0] = stag_lon[j, i]
    grid_corner_lon_local[iCell, 1] = stag_lon[j, i+1]
    grid_corner_lon_local[iCell, 2] = stag_lon[j+1, i+1]
    grid_corner_lon_local[iCell, 3] = stag_lon[j+1, i]
    grid_corner_lat_local[iCell, 0] = stag_lat[j, i]
    grid_corner_lat_local[iCell, 1] = stag_lat[j, i+1]
    grid_corner_lat_local[iCell, 2] = stag_lat[j+1, i+1]
    grid_corner_lat_local[iCell, 3] = stag_lat[j+1, i]

grid_corner_lon[:] = grid_corner_lon_local[:]
grid_corner_lat[:] = grid_corner_lat_local[:]

grid_imask[:] = 1  # For now, assume we don't want to mask anything out - but eventually may want to exclude certain cells from the input mesh during interpolation
grid_dims[:] = (nx * ny)


# plot some stuff
#import matplotlib.pyplot as plt
##plot a single point
#i=-1
#plt.plot(grid_center_lon[i], grid_center_lat[i], 'o')
#plt.plot(grid_corner_lon[i, 0], grid_corner_lat[i, 0], 'kx')
#plt.plot(grid_corner_lon[i, 1], grid_corner_lat[i, 1], 'bx')
#plt.plot(grid_corner_lon[i, 2], grid_corner_lat[i, 2], 'cx')
#plt.plot(grid_corner_lon[i, 3], grid_corner_lat[i, 3], 'gx')
#plt.show()
##plot all points
#plt.plot(grid_center_lon[:], grid_center_lat[:], 'bo')
#plt.plot(grid_corner_lon[:], grid_corner_lat[:], 'g.')
#plt.show()


fin.close()
fout.close()
