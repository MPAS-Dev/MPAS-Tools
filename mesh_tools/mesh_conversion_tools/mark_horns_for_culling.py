#!/usr/bin/env python
'''
This script identifies "horns" on a mesh (cells with two or fewer neighbors),
and marks them for culling.  In some cores/configurations, these weakly
connected cells can be dynamically inactive, and, therefore, undesirable to
keep in a mesh.

The method used will work on both planar and spherical meshes.
It adds the new masked cell to an existing 'cullCell' field if it exists,
otherwise it creates a new field.
'''

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import sys
import numpy as np
import netCDF4
from optparse import OptionParser
from datetime import datetime


print("== Gathering information.  (Invoke with --help for more details. All "
      "arguments are optional)\n")
parser = OptionParser()
parser.description = __doc__
parser.add_option(
    "-f",
    "--file",
    dest="inputFile",
    help="Name of file to be processed.",
    default="grid.nc",
    metavar="FILENAME")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

print("  File to be modified:  " + options.inputFile)


# Open file and get needed fields.
inputFile = netCDF4.Dataset(options.inputFile, 'r+')
nCells = len(inputFile.dimensions['nCells'])
cellsOnCell = inputFile.variables['cellsOnCell'][:]

# Add the horn cells to existing mask if it exists
if 'cullCell' in inputFile.variables:
    cullCell = inputFile.variables['cullCell'][:]
else:  # otherwise make a new mask initialized empty
    cullCell = np.zeros((nCells,))  # local variable

nHorns = 0
for i in range(nCells):
    # NOTE: Can change this threshold, if needed for a particular use case.
    if (cellsOnCell[i, :] > 0).sum() <= 2:
        cullCell[i] = 1
        nHorns += 1

# Write out the new field
if 'cullCell' in inputFile.variables:
    cullCellVar = inputFile.variables['cullCell']
else:
    cullCellVar = inputFile.createVariable('cullCell', 'i', ('nCells',))
cullCellVar[:] = cullCell


# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + \
    ": " + " ".join(sys.argv[:])
if hasattr(inputFile, 'history'):
    newhist = '\n'.join([thiscommand, getattr(inputFile, 'history')])
else:
    newhist = thiscommand
setattr(inputFile, 'history', newhist)

inputFile.close()

print('\n{} "horn" locations have been marked in the field cullCell.'.format(
    nHorns))
print("Remember to use MpasCellCuller.x to actually remove them!")
