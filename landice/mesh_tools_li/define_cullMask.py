#!/usr/bin/env python
# Script for adding a field named cullMask to an MPAS land ice grid for use with the MpasCellCuller tool that actually culls the unwanted cells.
# Matt Hoffman, February 28, 2013

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import sys
import numpy as np
from optparse import OptionParser
from netCDF4 import Dataset as NetCDFFile
from datetime import datetime


print("** Gathering information.")
parser = OptionParser()
parser.add_option("-f", "--file", dest="file", help="grid file to modify; default: landice_grid.nc", metavar="FILE")
parser.add_option("-m", "--method", dest="method", help="method to use for marking cells to cull.  Supported methods: 'noIce', 'numCells', 'distance', 'radius', 'edgeFraction'", metavar="METHOD")
parser.add_option("-n", "--numCells", dest="numCells", default=5, help="number of cells to keep beyond ice extent", metavar="NUM")
parser.add_option("-d", "--distance", dest="distance", default=50, help="distance (km) beyond ice extent to keep", metavar="DIST")
parser.add_option("-p", "--plot", dest="makePlot", help="Include to have the script generate a plot of the resulting mask, default=false", default=False, action="store_true")
options, args = parser.parse_args()

if not options.file:
        print("No grid filename provided. Using landice_grid.nc.")
        options.file = "landice_grid.nc"

if not options.method:
        sys.exit("ERROR: No method selected for choosing cells to mark for culling")
else:
        maskmethod = options.method

try:
  f = NetCDFFile(options.file,'r+')
except:
  sys.exit('Error loading grid file.')


xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
nCells = len(f.dimensions['nCells'])


# -- Get needed fields from the file --

cullCell = np.zeros((nCells, ), dtype=np.int8)  # Initialize to cull no cells


thicknessMissing = True
try:
  thickness = f.variables['thickness'][0,:]
  print('Using thickness field at time 0')
  thicknessMissing = False
except:
  print("The field 'thickness' is not available.  Some culling methods will not work.")


# =====  Various methods for defining the mask ====

# =========
#  only keep cells with ice
if maskmethod == 'noIce':
  print("Method: remove cells without ice")
  if thicknessMissing:
     sys.exit("Unable to perform 'numCells' method because thickness field was missing.")

  cullCell[thickness == 0.0] = 1

# =========
#  add a buffer of X cells around the ice
elif maskmethod == 'numCells':
  print("Method: remove cells beyond a certain number of cells from existing ice")

  if thicknessMissing:
     sys.exit("Unable to perform 'numCells' method because thickness field was missing.")

  buffersize=int(options.numCells)  # number of cells to expand
  print("Using a buffer of {} cells".format(buffersize))

  keepCellMask = np.copy(cullCell[:])
  keepCellMask[:] = 0
  cellsOnCell = f.variables['cellsOnCell'][:]
  nEdgesOnCell = f.variables['nEdgesOnCell'][:]

  # mark the cells with ice first
  keepCellMask[thickness > 0.0] = 1
  print('Num of cells with ice: {}'.format(sum(keepCellMask)))

  for i in range(buffersize):
    print('Starting buffer loop {}'.format(i+1))
    keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original
    ind = np.nonzero(keepCellMask == 0)[0]
    for i in range(len(ind)):
       iCell = ind[i]
       keepCellMaskNew[iCell] = keepCellMask[cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1].max() # if any neighbor has a value of 1, then 1 will get assigned to iCell.
    keepCellMask = np.copy(keepCellMaskNew)  # after we've looped over all cells assign the new mask to the variable we need (either for another loop around the domain or to write out)
    print('  Num of cells to keep: {}'.format(sum(keepCellMask)))

  # Now convert the keepCellMask to the cullMask
  cullCell[:] = np.absolute(keepCellMask[:]-1)  # Flip the mask for which ones to cull

# =========
#  remove cells beyond a certain distance of ice extent
elif maskmethod == 'distance':

  print("Method: remove cells beyond a certain distance from existing ice")

  if thicknessMissing:
     sys.exit("Unable to perform 'numCells' method because thickness field was missing.")

  dist=float(options.distance)
  print("Using a buffer distance of {} km".format(dist))
  dist = dist * 1000.0 # convert to m

  keepCellMask = np.copy(cullCell[:])
  keepCellMask[:] = 0
  cellsOnCell = f.variables['cellsOnCell'][:]
  nEdgesOnCell = f.variables['nEdgesOnCell'][:]
  xCell = f.variables['xCell'][:]
  yCell = f.variables['yCell'][:]

  # mark the cells with ice first
  keepCellMask[thickness > 0.0] = 1
  print('Num of cells with ice: {}'.format(sum(keepCellMask)))

  # find list of margin cells
  iceCells = np.nonzero(keepCellMask == 1)[0]
  marginMask = np.zeros((nCells, ), dtype=np.int8)
  for i in range(len(iceCells)):
     iCell = iceCells[i]
     for neighbor in cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1:  # the -1 converts from the fortran indexing in the variable to python indexing
         if thickness[neighbor] == 0.0:
            marginMask[iCell] = 1
            continue  # stop searching neighbors

  # loop over margin cells
  marginCells = np.nonzero(marginMask == 1)[0]
  for i in range(len(marginCells)):
      iCell = marginCells[i]
      # for each margin cell, find all cells within specified distance
      ind = np.nonzero(((xCell-xCell[iCell])**2 + (yCell-yCell[iCell])**2)**0.5 < dist)[0]
      keepCellMask[ind] = 1

  print('  Num of cells to keep:'.format(sum(keepCellMask)))

  # Now convert the keepCellMask to the cullMask
  cullCell[:] = np.absolute(keepCellMask[:]-1)  # Flip the mask for which ones to cull


# =========
#  cut out beyond some radius (good for the dome)
elif maskmethod == 'radius':
  print("Method: remove cells beyond a radius")
  ind = np.nonzero( (xCell[:]**2 + yCell[:]**2)**0.5 > 26000.0 )
  cullCell[ind] = 1

# =========
#  cut off some fraction of the height/width on all 4 sides - useful for cleaning up a mesh from periodic_general
elif maskmethod == 'edgeFraction':
  print("Method: remove a fraction from all 4 edges")
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


print('Num of cells to cull: {}'.format(sum(cullCell[:])))

# =========
# Try to add the new variable
if 'cullCell' not in f.variables:
  datatype = f.variables['indexToCellID'].dtype  # Get the datatype for integer
  f.createVariable('cullCell', datatype, ('nCells',))
f.variables['cullCell'][:] = cullCell

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
if hasattr(f, 'history'):
   newhist = '\n'.join([thiscommand, getattr(f, 'history')])
else:
   newhist = thiscommand
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
print("cullMask generation complete.")


