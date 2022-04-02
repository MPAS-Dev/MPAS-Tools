#!/usr/bin/env python
# Create a SCRIP file from a planar rectanfular mesh.
# See for details: http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_5_2_0rp1/ESMF_refdoc/node3.html#SECTION03024000000000000000

import sys
import netCDF4
import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
from pyproj import Transformer, transform, CRS

# ======== DEFINE PROJECTIONS =============
# Create empty dictionary to store projection definitions:
projections = dict()
# add more as needed:

# CISM's projection is as follows, with the vertical datum as EIGEN-GL04C geoid. 
# datum is actually EIGEN-GL04C but that is not an option in Proj.  Therefore using EGM08 which should be within ~1m everywhere (and 10-20 cm in most places)
# NOTE!!!!!!  egm08_25.gtx can be downloaded from:  http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx  and the path in the projection specification line should point to it!
#projections['gis-bamber'] = '+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +ellps=WGS84 +geoidgrids=./egm08_25.gtx'
projections['gis-bamber'] = '+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +ellps=WGS84' # This version ignores the vertical datum shift, which should be a very small error for horizontal-only positions

# GIMP projection: This is also polar stereographic but with different standard parallel and using the WGS84 ellipsoid.
projections['gis-gimp'] = '+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'

# BEDMAP2 projection
projections['ais-bedmap2'] = '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84' # Note: BEDMAP2 elevations use EIGEN-GL04C geoid

# Standard Lat/Long
projections['latlon'] = '+proj=longlat +ellps=WGS84'



print ("== Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser()
parser.description = "This script takes an MPAS grid file and generates a SCRIP grid file."
parser.add_option("-i", "--input", dest="inputFile", help="input grid file name used as input.", default="input.nc", metavar="FILENAME")
parser.add_option("-s", "--scrip", dest="scripFile", help="SCRIP grid file to output.", default="scrip.nc", metavar="FILENAME")
parser.add_option("-p", "--proj", dest="projection", help="projection used by the input data file. Valid options are:"  + str(projections.keys()), metavar="PROJ")
parser.add_option("-r", "--rank", dest="gridRank", help="desired rank of the output SCRIP grid data")
parser.add_option("--plot", dest="plot", action="store_true", help="if this flag is used, destination grid points are plotted")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

if not options.inputFile:
    sys.exit('Error: Data input grid file is required.  Specify with -c command line argument.')
if not options.scripFile:
    sys.exit('Error: SCRIP output grid file is required.  Specify with -s command line argument.')
if not options.projection:
    sys.exit('Error: data projection required with -p or --proj command line argument. Valid options are: ' + str(projections.keys()))
if not options.gridRank:
    sys.exit('Error: desired rank of SCRIP output grid data is required. Valid options are 1 (for unstructured grid) or 2')
print ('') # make a space in stdout before further output

# ===================================

fin = netCDF4.Dataset(options.inputFile, 'r')
fout = netCDF4.Dataset(options.scripFile, 'w')  # This will clobber existing files

# Get info from input file
x = fin.variables['x'][:]
y = fin.variables['y'][:]
nx = x.size
ny = y.size
dx = x[1] - x[0]
dy = y[1] - y[0]

# Write to output file
# Dimensions
fout.createDimension("grid_size", nx * ny)
fout.createDimension("grid_corners", 4 )

if int(options.gridRank) == 1:
    print('grid rank is 1')
    fout.createDimension("grid_rank", 1)
elif int(options.gridRank) == 2:
    print('grid rank is 2')
    fout.createDimension("grid_rank", 2)
else:
    raise ValueError(f'grid rank value is invalid: valid options are 1 or 2 but {options.gridRank} was given.')

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
print ('Building matrix version of x, y locations.')
xmatrix, ymatrix = np.meshgrid(x, y)
xc = np.append(x[:] - dx/2.0, x[-1] + dx/2.0 )  # get a copy of x that is on the staggered grid and includes both bounding edges
yc = np.append(y[:] - dy/2.0, y[-1] + dy/2.0 )  # get a copy of y that is on the staggered grid and includes both bounding edges
xcmatrix, ycmatrix = np.meshgrid(xc, yc)

# Unproject to lat/long for grid centers and grid corners
print ('Unprojecting.')

xmatrix_flat = xmatrix.flatten(order='C')  # Flatten using C indexing
ymatrix_flat = ymatrix.flatten(order='C')

# make a CRS (coordinate reference system) for projections from the Proj string:
crs_in = CRS.from_proj4(projections[options.projection])
crs_out = CRS.from_proj4(projections['latlon'])

# building a transformer
t = Transformer.from_crs(crs_in,crs_out)

# transform the original grid into the lat-lon grid
grid_center_lon[:], grid_center_lat[:] = t.transform(xmatrix_flat, ymatrix_flat)


# Now fill in the corners in the right locations
stag_lon, stag_lat = t.transform(xcmatrix, ycmatrix)
print ('Filling in corners of each cell.')
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

# set the grid dimension based on the grid rank
if int(options.gridRank) == 1:
    grid_dims[:] = (nx * ny)
elif int(options.gridRank) == 2:
    grid_dims[:] = [nx , ny]


if options.plot:
    print("plotting is on")
    # plot some stuff
    #plot a single point
    i=-1
    plt.figure(1)
    plt.plot(grid_center_lon[i], grid_center_lat[i], 'o')
    plt.plot(grid_corner_lon[i, 0], grid_corner_lat[i, 0], 'kx')
    plt.plot(grid_corner_lon[i, 1], grid_corner_lat[i, 1], 'bx')
    plt.plot(grid_corner_lon[i, 2], grid_corner_lat[i, 2], 'cx')
    plt.plot(grid_corner_lon[i, 3], grid_corner_lat[i, 3], 'gx')

    #plot all points
    plt.figure(2)
    plt.plot(grid_center_lon[:], grid_center_lat[:], 'bo')
    plt.plot(grid_corner_lon[:], grid_corner_lat[:], 'g.')
    plt.show()


fin.close()
fout.close()
