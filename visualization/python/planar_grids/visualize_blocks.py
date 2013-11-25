#!/usr/bin/python
# Basic script for visualizing block decompositionss.  Invoke with 'python plot_mpas_field.py --help for details about how to use.
# Matt Hoffman, June 14, 2012

import sys
import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
import netCDF4

print "** Gathering information."
parser = OptionParser()
parser.add_option("-g", "--gridfile", dest="gridfile", help="grid file to visualize; default: grid.nc", metavar="FILE")
parser.add_option("-b", "--blockfile", dest="blockfile", help="block file to visualize; default: graph.info.part.2", metavar="FILE")


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
  xCell = f.variables['xCell']
  yCell = f.variables['yCell']
except:
  sys.exit('xCell and/or yCell is/are missing.')

try:
  blocks = np.genfromtxt(options.blockfile)
except:
  sys.exit('block file could not be read properly.')

# Get the needed slice and determine the plot title.  Make some assumptions about how dimensions are arranged:
plottitle = 'block decomposition for grid file ' + options.gridfile + ' and graph file ' + options.blockfile


# MAKE THE PLOT
print '** Beginning to create plot.'
fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(xCell[:], yCell[:], 60, blocks, edgecolors='none') 
plt.colorbar()
plt.title( plottitle )

plt.draw()
plt.show()

f.close()

