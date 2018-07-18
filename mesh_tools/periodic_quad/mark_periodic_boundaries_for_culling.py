#!/usr/bin/env python
"""
This script will mark for culling the periodic rows/columns at the top, bottom,
left, and right of a mesh generated by periodic_quad.  It creates a new field
called cullCell and marks those locations with 1s.  This mesh can then be
run through mesh_conversion_tools/MpasCellCuller.x to remove those boundary
rows/columns, leaving a non-periodic planar mesh.

The logic to do this is to mark the max and min rows and
columns, so it is specific to a periodic_quad mesh.  If it were wished to
make this tool general, it could calculate cell-to-cell distances between
all neighboring cells and mark for culling any cell pairs that have
distances greater than, say, half of the range in x/y values of the entire mesh.
"""

import sys
import netCDF4
import numpy as np

from optparse import OptionParser


print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.description = "This script takes an MPAS grid file and marks the edge rows and columns for culling, e.g., to remove periodicity."
parser.add_option("-f", "--file", dest="inFile", help="MPAS grid file name used as input.", default="grid.nc", metavar="FILENAME")
parser.add_option("--keepx", dest="keepx", help="Retain periodicity in x-direction.", action='store_true', default=False)
parser.add_option("--keepy", dest="keepy", help="Retain periodicity in y-direction.", action='store_true', default=False)

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

print '' # make a space in stdout before further output


# ===============================================

fin = netCDF4.Dataset(options.inFile, 'r+')

# Get info from input file
xCell = fin.variables['xCell'][:]
yCell = fin.variables['yCell'][:]
nCells = len(fin.dimensions['nCells'])

if 'cullCell' in fin.variables:
  cullCell = fin.variables['cullCell']
else:
  cullCell = fin.createVariable('cullCell', fin.variables['xCell'].dtype, ('nCells',))

cullCell_local = np.zeros( (nCells,) )

# For a periodic quad, the upper and lower rows need to be marked
if not options.keepy:
    print 'y-periodic cells marked for culling'
    cullCell_local[np.nonzero(yCell == yCell.min())] = 1
    cullCell_local[np.nonzero(yCell == yCell.max())] = 1

if not options.keepy:
    print 'x-periodic cells marked for culling'
    cullCell_local[np.nonzero(xCell == xCell.min())] = 1
    cullCell_local[np.nonzero(xCell == xCell.max())] = 1

cullCell[:] = cullCell_local

fin.close()
print "Marking complete.  Don't forget to run MpasCellCuller!"
