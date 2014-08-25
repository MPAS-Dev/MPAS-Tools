#!/usr/bin/env python
import sys, os, glob, shutil, numpy, math

from netCDF4 import *
from netCDF4 import Dataset as NetCDFFile

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="Path to grid file", metavar="FILE")
parser.add_option("-s", "--scale", dest="scale", help="linear scale factor", default=1.0, metavar="FILE")
options, args = parser.parse_args()

if not options.filename:
	parser.error("A grid file is required.")

print "Applying scale factor of: ", options.scale

scale = float(options.scale)

grid = NetCDFFile(options.filename, 'a')

grid.variables['xCell'][:] = grid.variables['xCell'][:] * scale
grid.variables['yCell'][:] = grid.variables['yCell'][:] * scale
grid.variables['zCell'][:] = grid.variables['zCell'][:] * scale

grid.variables['xEdge'][:] = grid.variables['xEdge'][:] * scale
grid.variables['yEdge'][:] = grid.variables['yEdge'][:] * scale
grid.variables['zEdge'][:] = grid.variables['zEdge'][:] * scale

grid.variables['xVertex'][:] = grid.variables['xVertex'][:] * scale
grid.variables['yVertex'][:] = grid.variables['yVertex'][:] * scale
grid.variables['zVertex'][:] = grid.variables['zVertex'][:] * scale

grid.variables['dcEdge'][:] = grid.variables['dcEdge'][:] * scale
grid.variables['dvEdge'][:] = grid.variables['dvEdge'][:] * scale

grid.variables['areaCell'][:] = grid.variables['areaCell'][:] * scale**2
grid.variables['areaTriangle'][:] = grid.variables['areaTriangle'][:] * scale**2
grid.variables['kiteAreasOnVertex'][:] = grid.variables['kiteAreasOnVertex'][:] * scale**2

grid.close()
