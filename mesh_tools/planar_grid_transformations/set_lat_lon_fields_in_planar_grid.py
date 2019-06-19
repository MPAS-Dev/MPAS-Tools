#!/usr/bin/env python
'''
Take MPAS planar grid and populate the lat/lon fields based on a specified projection.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import netCDF4
import pyproj
from optparse import OptionParser
from datetime import datetime


# ======== DEFINE PROJECTIONS =============
# Create empty dictionary to store projection definitions:
projections = dict()
# add more as needed:

# CISM's projection is as follows, with the vertical datum as EIGEN-GL04C geoid.
# datum is actually EIGEN-GL04C but that is not an option in Proj.  Therefore using EGM08 which should be within ~1m everywhere (and 10-20 cm in most places)
# NOTE!!!!!!  egm08_25.gtx can be downloaded from:  http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx  and the path in the projection specification line should point to it!
#projections['gis-bamber'] = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids=./egm08_25.gtx')
projections['gis-bamber'] = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +ellps=WGS84')  # This version ignores the vertical datum shift, which should be a very small error for horizontal-only positions
projections['gis-bamber-shift'] = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=1300000.0 +y_0=3500000.0 +ellps=WGS84')  # This version ignores the vertical datum shift, which should be a very small error for horizontal-only positions

# GIMP projection: This is also polar stereographic but with different standard parallel and using the WGS84 ellipsoid.
projections['gis-gimp'] = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

# BEDMAP2 projection
projections['ais-bedmap2'] = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')  # Note: BEDMAP2 elevations use EIGEN-GL04C geoid

# Standard Lat/Long
projections['latlon'] = pyproj.Proj(proj='latlong', datum='WGS84')
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

print("Using {} projection, defined as: {}".format(options.projection, projections[options.projection].srs))

# get needed fields
f = netCDF4.Dataset(options.fileInName, 'r+')
xCell = f.variables['xCell']
yCell = f.variables['yCell']
xVertex = f.variables['xVertex']
yVertex = f.variables['yVertex']
xEdge = f.variables['xEdge']
yEdge = f.variables['yEdge']

latCell = f.variables['latCell']
lonCell = f.variables['lonCell']
latVertex = f.variables['latVertex']
lonVertex = f.variables['lonVertex']
latEdge = f.variables['latEdge']
lonEdge = f.variables['lonEdge']

print("Input file xCell min/max values:", xCell[:].min(), xCell[:].max())
print("Input file yCell min/max values:", yCell[:].min(), yCell[:].max())

# populate x,y fields
# MPAS uses lat/lon in radians, so have pyproj return fields in radians.
lonCell[:], latCell[:] = pyproj.transform(projections[options.projection], projections['latlon'], xCell[:], yCell[:], radians=True)
lonVertex[:], latVertex[:] = pyproj.transform(projections[options.projection], projections['latlon'], xVertex[:], yVertex[:], radians=True)
lonEdge[:], latEdge[:] = pyproj.transform(projections[options.projection], projections['latlon'], xEdge[:], yEdge[:], radians=True)

print("Calculated latCell min/max values (radians):", latCell[:].min(), latCell[:].max())
print("Calculated lonCell min/max values (radians):", lonCell[:].min(), lonCell[:].max())

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
if hasattr(f, 'history'):
   newhist = '\n'.join([thiscommand, getattr(f, 'history')])
else:
   newhist = thiscommand
setattr(f, 'history', newhist )


f.close()

print("Lat/lon calculations completed.  File has been written.")
