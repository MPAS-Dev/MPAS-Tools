#!/usr/bin/env python
'''
Matt Hoffman, 9/19/2022
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt

parser = OptionParser(description=__doc__)
parser.add_option("-n", dest="fileRegions", help="region file name.", metavar="FILENAME")
options, args = parser.parse_args()


f = Dataset(options.fileRegions, 'r')
regionCellMasks = f.variables['regionCellMasks'][:]
nRegions = len(f.dimensions['nRegions'])
nCells = len(f.dimensions['nCells'])

fout = Dataset("von_mises_calving_parameters.nc", 'w')
fout.createDimension('nCells', nCells)
fout.createDimension('Time', None)
grdVM = fout.createVariable('groundedVonMisesThresholdStress', 'd', ('Time', 'nCells',))
fltVM = fout.createVariable('floatingVonMisesThresholdStress', 'd', ('Time', 'nCells',))

values=[
        125.0,
        200.0,
        150.0,
        300.0, #300

        200.0, #225?
        350.0, # 400?
        400.0,
        125.0, # 130-400

        300.0, #400?
        300.0, #300-400
        300.0, # 300
        200.0, # 100?

        125.0,
        125.0,#?
        125.0,#120-400
        125.0,#125-300
        ]

grdVM[:]=100.0e3
fltVM[:]=100.0e3
for r in range(nRegions):
    mask = np.nonzero(regionCellMasks[:,r] == 1)[0]
    grdVM[0, mask] = values[r] * 1000.0
    fltVM[0, mask] = values[r] * 1000.0



fout.close()
f.close()
