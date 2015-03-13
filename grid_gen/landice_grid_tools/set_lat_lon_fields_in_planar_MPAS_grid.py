#!/usr/bin/env python
# Take MPAS planar land ice file and 1) shift its coordinates to a desired location and 2) populate the lat/lon fields

import sys
import netCDF4
import pyproj
#import numpy as np

from optparse import OptionParser


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

# Standard Lat/Long
projections['latlon'] = pyproj.Proj(proj='latlong', datum='WGS84')
# ===================================



print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.description = "This script 1) translates the coordinate system of the MPAS mesh specified with the -f flag so that the center of the domain is the center of the domain in the file specified by the -d flag; 2) populates all of the MPAS lat and lon fields based on the projection specified by the -p option."
parser.add_option("-f", "--file", dest="fileInName", help="MPAS land ice file name.", default="landice_grid.nc", metavar="FILENAME")
parser.add_option("-d", "--datafile", dest="dataFileName", help="data file name from which to get domain information.  Required.  Uses x1 and y1 fields.", metavar="FILENAME")
parser.add_option("-p", "--proj", dest="projection", help="projection used for the data. Valid options are:"  + str(projections.keys()), metavar="PROJ")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

if not options.projection:
    sys.exit('Error: data projection required with -p or --proj command line argument. Valid options are: ' + str(projections.keys()))
if not options.dataFileName:
    sys.exit('Error: datafile required with -d or --datafile command line argument.')

if not options.fileInName:
    print "No output filename specified, so using 'landice_grid.nc'."
    options.fileInName = 'landice_grid.nc'
print '' # make a space in stdout before further output


# =================================================

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

fd = netCDF4.Dataset(options.dataFileName, 'r')
x1 = fd.variables['x1'][:]
y1 = fd.variables['y1'][:]


# step 1, perform a shift

mpasXcenter = (xCell[:].min() + xCell[:].max()) * 0.5
mpasYcenter = (yCell[:].min() + yCell[:].max()) * 0.5


dataXcenter = (x1[:].min() + x1[:].max()) * 0.5
dataYcenter = (y1[:].min() + y1[:].max()) * 0.5


xOffset = dataXcenter - mpasXcenter
yOffset = dataYcenter - mpasYcenter

xCell[:] += xOffset
yCell[:] += yOffset
xVertex[:] += xOffset
yVertex[:] += yOffset
xEdge[:] += xOffset
yEdge[:] += yOffset


# step 2, populate x,y fields
# MPAS uses lat/lon in radians, so have pyproj return fields in radians.
lonCell[:], latCell[:] = pyproj.transform(projections[options.projection], projections['latlon'], xCell[:], yCell[:], radians=True)
lonVertex[:], latVertex[:] = pyproj.transform(projections[options.projection], projections['latlon'], xVertex[:], yVertex[:], radians=True)
lonEdge[:], latEdge[:] = pyproj.transform(projections[options.projection], projections['latlon'], xEdge[:], yEdge[:], radians=True)


f.close()
fd.close()
