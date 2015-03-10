#!/usr/bin/env python
# Mark periodic rows/columns for deletion

import sys
import netCDF4
import numpy as np

from optparse import OptionParser


print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.description = "This script takes an MPAS grid file and marks the edge rows and columns for culling, e.g., to remove periodicity."
parser.add_option("-f", "--file", dest="inFile", help="MPAS grid file name used as input.", default="grid.nc", metavar="FILENAME")
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

# For a periodic hex, the upper and lower rows need to be marked
cullCell_local[np.nonzero(yCell == yCell.min())] = 1
cullCell_local[np.nonzero(yCell == yCell.max())] = 1

# For a periodidic hex the leftmost and rightmost *TWO* columns need to be marked
unique_Xs=np.array(sorted(list(set(xCell[:]))))
cullCell_local[np.nonzero(xCell == unique_Xs[0])] = 1
cullCell_local[np.nonzero(xCell == unique_Xs[1])] = 1
cullCell_local[np.nonzero(xCell == unique_Xs[-1])] = 1
cullCell_local[np.nonzero(xCell == unique_Xs[-2])] = 1

cullCell[:] = cullCell_local

fin.close()

