#!/usr/bin/python
# Script for adding a field named cullMask to an MPAS land ice grid for use with the crop_landice_grid.F tool that actually culls the unwanted cells.
# Matt Hoffman, February 28, 2013

import sys
import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
try:
  from Scientific.IO.NetCDF import NetCDFFile
  netCDF_module = 'Scientific.IO.NetCDF'
except ImportError:
  try:
    from netCDF4 import Dataset as NetCDFFile
    netCDF_module = 'netCDF4'
  except ImportError:
      print 'Unable to import any of the following python modules:'
      print '  Scientific.IO.NetCDF \n  netcdf4 '
      print 'One of them must be installed.'
      raise ImportError('No netCDF module found')


print "** Gathering information."
parser = OptionParser()
parser.add_option("-f", "--file", dest="file", help="grid file to modify; default: land_ice_grid.nc", metavar="FILE")


options, args = parser.parse_args()

if not options.file:
	print "No grid filename provided. Using land_ice_grid.nc."
        options.file = "land_ice_grid.nc"

try: 
  f = NetCDFFile(options.file,'r+')
except:
  sys.exit('Error loading grid file.')
  

# get the right representation of the datatype
if netCDF_module == 'Scientific.IO.NetCDF':
    datatype = 'i'
else:
    datatype = f.variables['indexToCellID'].dtype  # Get the datatype for double precision float

# Try to add the new variable
if 'keepCellMask' not in f.variables:
  f.createVariable('keepCellMask', datatype, ('nCells',))

keepCellMask = f.variables['keepCellMask'][:]
thickness = f.variables['thickness'][:]

keepCellMask[:] = 0

# =====  Various methods for defining the mask ====
maskmeth = 2

# 1. only keep cells with ice
if maskmeth == 1:
  keepCellMask[thickness[0,:] > 0.0] = 1

# 2. add a buffer of X cells around the ice 
elif maskmeth == 2:

  buffersize=4  # number of cells to expand

  keepCellMask[thickness[0,:] > 0.0] = 1  # start with the mask being where ice is
  print 'Num of cells to keep:', sum(keepCellMask)
  cellsOnCell = f.variables['cellsOnCell']
  for i in range(buffersize):
    print 'Starting buffer loop ', i+1
    keepCellMaskNew = keepCellMask[:]+0 # make a copy to edit  --- the +0 makes an actual copy rather than pointing to the memory
    for iCell in range(len(keepCellMask)):
      for neighbor in cellsOnCell[iCell,:]-1:  # the -1 converts from the fortran indexing in the variable to python indexing
        if keepCellMask[neighbor] == 1:  # if any neighbors are already being kept on the old mask then keep this cell too.  This will get all (or almost all) of the ice cells, but they have already been assigned so they don't matter.  What we care about is getting cells that are not ice that have on or more ice neighbors.
           keepCellMaskNew[iCell] = 1
    keepCellMask = keepCellMaskNew[:]+0  # after we've looped over all cells assign the new mask to the variable we need (either for another loop around the domain or to write out)
    print 'Num of cells to keep:', sum(keepCellMask)
fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(f.variables['xCell'][:], f.variables['yCell'][:], 50, keepCellMask, edgecolors='none') #, vmin=minval, vmax=maxval)
plt.colorbar()
plt.draw()
plt.show()

f.variables['keepCellMask'][:] = keepCellMask

#f.close()



