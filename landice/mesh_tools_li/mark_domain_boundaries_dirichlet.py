#!/usr/bin/env python
'''
This script marks all of the boundary cells in a domain as Dirichlet velocity boundaries.
'''

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import netCDF4
import numpy as np
from optparse import OptionParser
from datetime import datetime
import sys

print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n")
parser = OptionParser()
parser.description = __doc__
parser.add_option("-f", "--file", dest="inputFile", help="name of file to be modified.", default="landice_grid.nc", metavar="FILENAME")
parser.add_option("-t", "--time", dest="time", help="time level to modify", default=0, type="int",  metavar="FILENAME")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()


print("  Input file: {}".format(options.inputFile))
print("  Time level: {}".format(options.time))

f=netCDF4.Dataset(options.inputFile, 'r+')
nCells = len(f.dimensions['nCells'])
mask = f.variables['dirichletVelocityMask'][options.time, :, :]
cONc = f.variables['cellsOnCell'][:]
nEdgesOnCell = f.variables['nEdgesOnCell'][:]

mask[:] = 0
for i in range(nCells):
   nE = nEdgesOnCell[i]
   if min(cONc[i, :nE]) == 0:
      mask[i,:] = 1
f.variables['dirichletVelocityMask'][options.time, :, :] = mask[:]


# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
if hasattr(f, 'history'):
   newhist = '\n'.join([thiscommand, getattr(f, 'history')])
else:
   newhist = thiscommand
setattr(f, 'history', newhist )

f.close()

print('\nMarking boundary cells completed.')
