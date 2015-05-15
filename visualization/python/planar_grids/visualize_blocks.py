#!/usr/bin/env python
# Basic script for visualizing block decompositionss.  Invoke with 'python plot_mpas_field.py --help for details about how to use.
# Matt Hoffman, June 14, 2012

import sys
import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
import netCDF4
import matplotlib.cm as cm
import matplotlib

print "** Gathering information."
parser = OptionParser()
parser.add_option("-g", "--gridfile", dest="gridfile", help="grid file to visualize; default: grid.nc", metavar="FILE")
parser.add_option("-b", "--blockfile", dest="blockfile", help="block file to visualize; default: graph.info.part.2", metavar="FILE")
parser.add_option("-c", "--colormap", dest="cmap", help="color map to use: 'jet' or 'random'; default: 'random'", metavar="COLORMAP")
parser.add_option("-m", "--mark", dest="mark", help="block number to mark with black stars", metavar="NUMBER")


options, args = parser.parse_args()

if not options.gridfile:
	print "No grid filename provided. Using grid.nc."
        options.gridfile = "grid.nc"

if not options.blockfile:
	print "No block filename provided. Using graph.info.part.2."
        options.blockfile = "graph.info.part.2"

try:
  f = netCDF4.Dataset(options.gridfile,'r')
except:
  sys.exit('Error loading grid file.')


# Get grid stuff
#times = f.variables['xtime']  # Not needed unless trying to print actual time stamp
try:
  xCell = f.variables['xCell'][:]
  yCell = f.variables['yCell'][:]
  nCells = len(f.dimensions['nCells'])
except:
  sys.exit('xCell and/or yCell is/are missing.')

try:
  blocks = np.genfromtxt(options.blockfile, dtype='int')
except:
  sys.exit('block file could not be read properly.')

if nCells != len(blocks):
    sys.exit('Error: Number of lines in block file does not equal nCells in the grid file!')

counts = np.bincount(blocks)
print '  Min number of cells per block: ', counts.min()
print '  Mean number of cells per block:', counts.mean()
print '  Max number of cells per block: ', counts.max()

# MAKE THE PLOT

if options.cmap == 'jet':
   cmap = plt.get_cmap('jet')
else:
   cmap = matplotlib.colors.ListedColormap ( np.random.rand (blocks.max(), 3))

print '** Beginning to create plot.'
plottitle = 'block decomposition for grid file ' + options.gridfile + '\nand graph file ' + options.blockfile
fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(xCell[:], yCell[:], c=blocks, s=12, edgecolors='none' , cmap=cmap)
plt.colorbar()
plt.title( plottitle )

if options.mark:
   ind = np.where(blocks == int(options.mark))
   plt.plot(xCell[ind], yCell[ind], 'k*')

plt.draw()
plt.show()

f.close()

