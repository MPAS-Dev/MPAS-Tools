#!/usr/bin/env python
# Script for adding a field named cullMask to an MPAS land ice grid for use with the crop_landice_grid.F tool that actually culls the unwanted cells.
# Matt Hoffman, February 28, 2013

import sys
import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
from netCDF4 import Dataset as NetCDFFile


print "** Gathering information."
parser = OptionParser()
parser.add_option("-f", "--file", dest="file", help="grid file to modify; default: landice_grid.nc", metavar="FILE")


options, args = parser.parse_args()

if not options.file:
	print "No grid filename provided. Using landice_grid.nc."
        options.file = "landice_grid.nc"

try: 
  f = NetCDFFile(options.file,'r+')
except:
  sys.exit('Error loading grid file.')
  

datatype = f.variables['indexToCellID'].dtype  # Get the datatype for integer

# Try to add the new variable
if 'cullCell' not in f.variables:
  f.createVariable('cullCell', datatype, ('nCells',))

# -- Get needed fields from the file --
cullCell = f.variables['cullCell']
cullCell[:] = 0  # Initialize to cull no cells

xCell = f.variables['xCell']
yCell = f.variables['yCell']

try:
  thickness = f.variables['thickness']
except:
  print "The field 'thickness' is not available.  Some culling methods will not work."


# =====  Various methods for defining the mask ====
maskmeth = 2

# 1. only keep cells with ice
if maskmeth == 1:
  cullCell[thickness[0,:] == 0.0] = 1

# 2. add a buffer of X cells around the ice 
elif maskmeth == 2:

  buffersize=2  # number of cells to expand

  keepCellMask = np.copy(cullCell[:])   # to use existing code, do this.  not very efficient but fine for now.
  keepCellMask[:] = 0
  keepCellMask[thickness[0,:] > 0.0] = 1  # start with the mask being where ice is
  print 'Num of cells to keep:', sum(keepCellMask)
  cellsOnCell = f.variables['cellsOnCell']
  for i in range(buffersize):
    print 'Starting buffer loop ', i+1
    keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original
    for iCell in range(len(keepCellMask)):
      if keepCellMask[iCell] == 0:  # don't bother to check cells we are already keeping
        for neighbor in cellsOnCell[iCell,:]-1:  # the -1 converts from the fortran indexing in the variable to python indexing
          if neighbor>-1 and keepCellMask[neighbor] == 1:  # if any neighbors are already being kept on the old mask then keep this cell too.  This will get all (or almost all) of the ice cells, but they have already been assigned so they don't matter.  What we care about is getting cells that are not ice that have one or more ice neighbors.
             keepCellMaskNew[iCell] = 1
    keepCellMask = np.copy(keepCellMaskNew)  # after we've looped over all cells assign the new mask to the variable we need (either for another loop around the domain or to write out)
    print 'Num of cells to keep:', sum(keepCellMask)
  # Now convert the keepCellMask to the cullMask
  cullCell[:] = np.absolute(keepCellMask[:]-1)  # Flip the mask for which ones to cull
  print 'Num of cells to cull:', sum(cullCell[:])


# 3. cut out beyond some radius (good for the dome)
elif maskmeth == 3:
  ind = np.nonzero( (xCell[:]**2 + yCell[:]**2)**0.5 > 28000.0 )
  cullCell[ind] = 1

###########
## SHIFT VALUES to 0
#xshift = -1 * xCell[:].min()
#yshift = -1 * yCell[:].min()
#xCell[:] = xCell[:] + xshift
#yCell[:] = yCell[:] + yshift
#f.variables['xVertex'][:] = f.variables['xVertex'][:] + xshift
#f.variables['yVertex'][:] = f.variables['yVertex'][:] + yshift
#f.variables['xEdge'][:] = f.variables['xEdge'][:] + xshift
#f.variables['yEdge'][:] = f.variables['yEdge'][:] + yshift

fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(xCell[:], yCell[:], 50, cullCell[:], edgecolors='none') #, vmin=minval, vmax=maxval)
plt.colorbar()
plt.draw()
plt.show()


f.close()



