#!/usr/bin/env python
'''
Translate planar MPAS grid by one of three methods
'''

import sys
import netCDF4
#import numpy as np
from optparse import OptionParser


print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.description = ("This script translates the coordinate system of the planar MPAS mesh specified with the -f flag. "
                      "There are 3 possible methods to choose from:"
                      "1) shift the origin to the center of the domain"
                      "2) arbirary shift in x and/or y"
                      "3) shift to the center of the domain described in a separate file")
parser.add_option("-f", "--file", dest="fileInName", help="MPAS planar grid file name.", default="grid.nc", metavar="FILENAME")
parser.add_option("-d", "--datafile", dest="dataFileName", help="data file name to which to match the domain center of.  Uses xCell,yCell or, if those fields do not exist, will secondly try x1,y1 fields.", metavar="FILENAME")
parser.add_option("-x", dest="xshift", help="user-specified shift in the x-direction.", type="float", default=0.0, metavar="SHIFT_VALUE")
parser.add_option("-y", dest="yshift", help="user-specified shift in the y-direction.", type="float", default=0.0, metavar="SHIFT_VALUE")
parser.add_option("-c", dest="center", help="shift so origin is at center of domain", action="store_true", default=False)
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

print "Attempting to translate coordinates in file: %s"%options.fileInName


if options.dataFileName and (options.xshift or options.yshift):
  sys.exit("Error: Specifying a datafile AND one or both of x/y shift is invalid.  Please select one of those methods only.")

if options.center and (options.xshift or options.yshift):
  sys.exit("Error: Specifying a shift to center AND one or both of x/y shift is invalid.  Please select one of those methods only.")

if options.dataFileName and options.center:
  sys.exit("Error: Specifying a datafile AND a shift to center is invalid.  Please select one of those methods only.")

if not options.center and not options.xshift and not options.yshift and not options.dataFileName:
  sys.exit("Error: No translation method was specified.  Please select one.  Run with -h for more information.")

if options.dataFileName:
  method = 'file'
  print "  Translating coordinates in %s so the domain center matches the domain center in %s."%(options.fileInName, options.dataFileName)

if options.xshift or options.yshift:
  method = 'xy'
  print "  Translating coordinates in %s by user-specified values.  X-shift=%f; Y-shift=%f"%(options.fileInName, options.xshift, options.yshift)

if options.center:        
  method = 'center'
  print "  Translating coordinates in %s so the origin is the center of the domain."

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

mpasXcenter = (xCell[:].min() + xCell[:].max()) * 0.5
mpasYcenter = (yCell[:].min() + yCell[:].max()) * 0.5


if method == 'file':
   fd = netCDF4.Dataset(options.dataFileName, 'r')
   try:
     x = fd.variables['xCell'][:]
     y = fd.variables['yCell'][:]
   except:
     try:
       x = fd.variables['x1'][:]
       y = fd.variables['y1'][:]
     except:
       sys.exit('Error: data file specified has neither xCell/yCell nor x1/y1 fields.')
   dataXcenter = (x[:].min() + x[:].max()) * 0.5
   dataYcenter = (y[:].min() + y[:].max()) * 0.5
   fd.close()
   xOffset = dataXcenter - mpasXcenter
   yOffset = dataYcenter - mpasYcenter
elif method == 'xy':
   xOffset = options.xshift
   yOffset = options.yshift
elif method == 'center':
   xOffset = -1.0 * mpasXcenter
   yOffset = -1.0 * mpasYcenter

# perform the shift
xCell[:] += xOffset
yCell[:] += yOffset
xVertex[:] += xOffset
yVertex[:] += yOffset
xEdge[:] += xOffset
yEdge[:] += yOffset

# Update history attribute of netCDF file
if hasattr(f, 'history'):
   newhist = '\n'.join([getattr(f, 'history'), ' '.join(sys.argv[:]) ] )
else:
   newhist = sys.argv[:]
setattr(f, 'history', newhist )

f.close()

print "Translation completed."
