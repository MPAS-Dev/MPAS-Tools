#!/usr/bin/env python
# Script for adding a field named cullMask to an MPAS land ice grid for use with the MpasCellCuller tool that actually culls the unwanted cells.
# Matt Hoffman, February 28, 2013

import sys
import numpy as np
from optparse import OptionParser
from netCDF4 import Dataset as NetCDFFile


print "** Gathering information."
parser = OptionParser()
parser.add_option("-f", "--file", dest="file", help="grid file to modify; default: landice_grid.nc", metavar="FILE")
parser.add_option("-m", "--method", dest="method", help="method to use for marking cells to cull.  See code for available options, most of which will need tweaking for particular applications", metavar="METHOD")
parser.add_option("-p", "--plot", dest="makePlot", help="Include to have the script generate a plot of the resulting mask, default=false", default=False, action="store_true")
options, args = parser.parse_args()

if not options.file:
	print "No grid filename provided. Using landice_grid.nc."
        options.file = "landice_grid.nc"

if not options.method:
        sys.exit("ERROR: No method selected for choosing cells to mark for culling")
else:
        maskmethod = int(options.method)

try:
  f = NetCDFFile(options.file,'r+')
except:
  sys.exit('Error loading grid file.')


xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
nCells = len(f.dimensions['nCells'])


# -- Get needed fields from the file --

cullCell = np.zeros((nCells, ), dtype=np.int8)  # Initialize to cull no cells


try:
  thickness = f.variables['thickness'][0,:]
  print 'Using thickness field at time 0'
except:
  print "The field 'thickness' is not available.  Some culling methods will not work."


# =====  Various methods for defining the mask ====

# =========
# 1. only keep cells with ice
if maskmethod == 1:
  print "Method: remove cells without ice"
  cullCell[thickness == 0.0] = 1

# =========
# 2. add a buffer of X cells around the ice 
elif maskmethod == 2:
  print "Method: remove cells beyond a certain number of cells from existing ice"

  buffersize=25  # number of cells to expand

  keepCellMask = np.copy(cullCell[:])
  keepCellMask[:] = 0
  cellsOnCell = f.variables['cellsOnCell'][:]
  nEdgesOnCell = f.variables['nEdgesOnCell'][:]

  # mark the cells with ice first
  keepCellMask[thickness > 0.0] = 1
  print 'Num of cells with ice:', sum(keepCellMask)

  for i in range(buffersize):
    print 'Starting buffer loop ', i+1
    keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original
    for iCell in range(len(keepCellMask)):
      if keepCellMask[iCell] == 0:  # don't bother to check cells we are already keeping
        for neighbor in cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1:  # the -1 converts from the fortran indexing in the variable to python indexing
          if neighbor >= 0 and keepCellMask[neighbor] == 1:  # if any neighbors are already being kept on the old mask then keep this cell too.  This will get a lot of the ice cells, but they have already been assigned so they don't matter.  What we care about is getting cells that are not ice that have one or more ice neighbors.
             keepCellMaskNew[iCell] = 1
    keepCellMask = np.copy(keepCellMaskNew)  # after we've looped over all cells assign the new mask to the variable we need (either for another loop around the domain or to write out)
    print '  Num of cells to keep:', sum(keepCellMask)

  # Now convert the keepCellMask to the cullMask
  cullCell[:] = np.absolute(keepCellMask[:]-1)  # Flip the mask for which ones to cull
  del keepCellMask, keepCellMaskNew

# =========
# 3. cut out beyond some radius (good for the dome)
elif maskmethod == 3:
  print "Method: remove cells beyond a radius"
  ind = np.nonzero( (xCell[:]**2 + yCell[:]**2)**0.5 > 28000.0 )
  cullCell[ind] = 1

# =========
# 4. cut off some fraction of the height/width on all 4 sides - useful for cleaning up a mesh from periodic_general
elif maskmethod == 4:
  print "Method: remove a fraction from all 4 edges"
  frac=0.025

  cullCell[:] = 0
  width = xCell.max()-xCell.min()
  height = yCell.max()-yCell.min()
  ind = np.nonzero( xCell[:] < (xCell.min() + width*frac) )
  cullCell[ind] = 1
  ind = np.nonzero( xCell[:] > (xCell.max() - width*frac) )
  cullCell[ind] = 1
  ind = np.nonzero( yCell[:] < (yCell.min() + height*frac) )
  cullCell[ind] = 1
  ind = np.nonzero( yCell[:] > (yCell.max() - height*frac) )
  cullCell[ind] = 1

# =========
else:
  sys.exit("ERROR: no valid culling method selected.")
# =========


print 'Num of cells to cull:', sum(cullCell[:])

# =========
# Try to add the new variable
if 'cullCell' not in f.variables:
  datatype = f.variables['indexToCellID'].dtype  # Get the datatype for integer
  f.createVariable('cullCell', datatype, ('nCells',))
f.variables['cullCell'][:] = cullCell

# Update history attribute of netCDF file
if hasattr(f, 'history'):
   newhist = '\n'.join([getattr(f, 'history'), ' '.join(sys.argv[:]) ] )
else:
   newhist = sys.argv[:]
setattr(f, 'history', newhist )


# --- make a plot only if requested ---
if options.makePlot:
  import matplotlib.pyplot as plt
  fig = plt.figure(1, facecolor='w')
  ax = fig.add_subplot(111, aspect='equal')
  plt.scatter(xCell[:], yCell[:], 50, cullCell[:], edgecolors='none') #, vmin=minval, vmax=maxval)
  plt.colorbar()
  plt.draw()
  plt.show()

f.close()
print "cullMask generation complete."


