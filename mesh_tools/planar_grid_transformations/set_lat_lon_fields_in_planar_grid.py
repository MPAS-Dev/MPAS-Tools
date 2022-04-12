#!/usr/bin/env python
'''
Take MPAS planar grid and populate the lat/lon fields based on a specified projection.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
import netCDF4
from pyproj import Transformer, transform, CRS
from optparse import OptionParser
from datetime import datetime


# ======== DEFINE PROJECTIONS =============
# Create empty dictionary to store projection definitions:
projections = dict()
# add more as needed:

# CISM's projection is as follows, with the vertical datum as EIGEN-GL04C geoid.
# datum is actually EIGEN-GL04C but that is not an option in Proj.  Therefore using EGM08 which should be within ~1m everywhere (and 10-20 cm in most places)
# NOTE!!!!!!  egm08_25.gtx can be downloaded from:  http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx  and the path in the projection specification line should point to it!
#projections['gis-bamber'] = '+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids=./egm08_25.gtx'
projections['gis-bamber'] = '+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +ellps=WGS84'  # This version ignores the vertical datum shift, which should be a very small error for horizontal-only positions
projections['gis-bamber-shift'] = '+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=1300000.0 +y_0=3500000.0 +ellps=WGS84'  # This version ignores the vertical datum shift, which should be a very small error for horizontal-only positions

# GIMP projection: This is also polar stereographic but with different standard parallel and using the WGS84 ellipsoid.
projections['gis-gimp'] = '+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'

# BEDMAP2 projection
projections['ais-bedmap2'] = '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'  # Note: BEDMAP2 elevations use EIGEN-GL04C geoid

# Standard Lat/Long
projections['latlon'] = '+proj=longlat +ellps=WGS84'
# ===================================



print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser()
parser.description = "This script populates the MPAS lat and lon fields based on the projection specified by the -p option."
parser.add_option("-f", "--file", dest="fileInName", help="MPAS land ice file name.", default="landice_grid.nc", metavar="FILENAME")
parser.add_option("-p", "--proj", dest="projection", help="projection used for the data. Valid options are: \n"  + str(list(projections.keys())), metavar="PROJ")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

if not options.projection:
    sys.exit('Error: data projection required with -p or --proj command line argument. Valid options are: ' + str(list(projections.keys())))

if not options.fileInName:
    print("No filename specified, so using 'landice_grid.nc'.")
    options.fileInName = 'landice_grid.nc'
print('') # make a space in stdout before further output


# =================================================

print("Using {} projection, defined as: {}".format(options.projection, projections[options.projection]))

# get needed fields
f = netCDF4.Dataset(options.fileInName, 'r+')
xCell = f.variables['xCell']
yCell = f.variables['yCell']
xVertex = f.variables['xVertex']
yVertex = f.variables['yVertex']
xEdge = f.variables['xEdge']
yEdge = f.variables['yEdge']

latCellVar = f.variables['latCell']
lonCellVar = f.variables['lonCell']
latVertexVar = f.variables['latVertex']
lonVertexVar = f.variables['lonVertex']
latEdgeVar = f.variables['latEdge']
lonEdgeVar = f.variables['lonEdge']

print("Input file xCell min/max values:", xCell[:].min(), xCell[:].max())
print("Input file yCell min/max values:", yCell[:].min(), yCell[:].max())

# make a CRS (coordinate reference system) for projections from Proj string:
crs_in = CRS.from_proj4(projections[options.projection])
crs_out = CRS.from_proj4(projections['latlon'])

# define a transformer
t = Transformer.from_crs(crs_in,crs_out)

# populate x,y fields
# MPAS uses lat/lon in radians, so have pyproj return fields in radians.
lonCell, latCell = t.transform(xCell[:], yCell[:], radians=True)
lonVertex, latVertex = t.transform(xVertex[:], yVertex[:], radians=True)
lonEdge, latEdge = t.transform(xEdge[:], yEdge[:], radians=True)

# change the longitude convention to use positive values [0 2pi]
lonCell = np.mod(lonCell, 2.0*np.pi)
lonVertex = np.mod(lonVertex, 2.0*np.pi)
lonEdge = np.mod(lonEdge, 2.0*np.pi)

print("Calculated latCell min/max values (radians):", latCell.min(), latCell.max())
print("Calculated lonCell min/max values (radians):", lonCell.min(), lonCell.max())

latCellVar[:] = latCell
lonCellVar[:] = lonCell
latVertexVar[:] = latVertex
lonVertexVar[:] = lonVertex
latEdgeVar[:] = latEdge
lonEdgeVar[:] = lonEdge

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
if hasattr(f, 'history'):
   newhist = '\n'.join([thiscommand, getattr(f, 'history')])
else:
   newhist = thiscommand
setattr(f, 'history', newhist )


f.close()

print("Lat/lon calculations completed.  File has been written.")
