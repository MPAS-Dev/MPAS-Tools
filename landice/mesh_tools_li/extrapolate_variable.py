#!/usr/bin/env python
"""
Script to extrapolate arbitrary variable

Created on Mon Feb1 2021

@author: Matt Hoffman, Trevor Hillebrand
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import sys
from datetime import datetime


parser = OptionParser(description='Convert data from exodus file to MPAS mesh. WARNING: Change the SEACAS library dir to your own path! A simple usage example: conversion_exodus_init_to_mpasli_mesh.py -e antarctica.exo -o target.nc -a ./mpas_cellID.ascii -v beta')

parser.add_option("-f", "--file", dest="nc_file", help="the mpas file to write to")
parser.add_option("-v", "--variable", dest="var_name", help="the mpas variable(s) you want to extrapolate")
parser.add_option("-m", "--method", dest="extrapolation", default='min', help="idw or min method of extrapolation")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

dataset = Dataset(options.nc_file, 'r+')
var_name = options.var_name
varValue = dataset.variables[var_name][0, :]
extrapolation = options.extrapolation
# Extrapolation
nCells = len(dataset.dimensions['nCells'])
if 'thickness' in dataset.variables.keys():
    thickness = dataset.variables['thickness'][0,:]
    bed = dataset.variables["bedTopography"][0,:]
cellsOnCell = dataset.variables['cellsOnCell'][:]
nEdgesOnCell = dataset.variables['nEdgesOnCell'][:]
xCell = dataset.variables["yCell"][:]
yCell = dataset.variables["xCell"][:]

# Define region of good data to extrapolate from.  Different methods for different variables
if var_name in ["effectivePressure", "beta", "muFriction"]:
    groundedMask = (thickness > (-1028.0 / 910.0 * bed))
    keepCellMask = np.copy(groundedMask)
    extrapolation == "min"

    for iCell in range(nCells):
        for n in range(nEdgesOnCell[iCell]):
            jCell = cellsOnCell[iCell, n] - 1
            if (groundedMask[jCell] == 1):
                keepCellMask[iCell] = 1
                continue
    keepCellMask *= (varValue > 0)  # ensure zero muFriction does not get extrapolated
elif var_name in ["floatingBasalMassBal"]:
    floatingMask = (thickness <= (-1028.0 / 910.0 * bed))
    keepCellMask = floatingMask * (varValue != 0.0)
    extrapolation == "idw"
else:
    keepCellMask = (thickness > 0.0)

keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that will be used later
keepCellMaskOrig = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original

# recursive extrapolation steps:
# 1) find cell A with mask = 0
# 2) find how many surrounding cells have nonzero mask, and their indices (this will give us the cells from upstream)
# 3) use the values for nonzero mask upstream cells to extrapolate the value for cell A
# 4) change the mask for A from 0 to 1
# 5) Update mask
# 6) go to step 1)

print("\nStart {} extrapolation using {} method".format(var_name, extrapolation))
while np.count_nonzero(keepCellMask) != nCells:

    keepCellMask = np.copy(keepCellMaskNew)
    searchCells = np.where(keepCellMask==0)[0]
    varValueOld = np.copy(varValue)

    for iCell in searchCells:
        neighborcellID = cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1
        neighborcellID = neighborcellID[neighborcellID>=0] # Important: ignore the phantom "neighbors" that are off the edge of the mesh (0 values in cellsOnCell)

        mask_for_idx = keepCellMask[neighborcellID] # active cell mask
        mask_nonzero_idx, = np.nonzero(mask_for_idx)

        nonzero_id = neighborcellID[mask_nonzero_idx] # id for nonzero beta cells
        nonzero_num = np.count_nonzero(mask_for_idx)

        assert len(nonzero_id) == nonzero_num

        if nonzero_num > 0:
            x_i = xCell[iCell]
            y_i = yCell[iCell]
            x_adj = xCell[nonzero_id]
            y_adj = yCell[nonzero_id]
            var_adj = varValueOld[nonzero_id]
            if extrapolation == 'idw':
                ds = np.sqrt((x_i-x_adj)**2+(y_i-y_adj)**2)
                assert np.count_nonzero(ds)==len(ds)
                var_interp = 1.0/sum(1.0/ds)*sum(1.0/ds*var_adj)
                varValue[iCell] = var_interp
            elif extrapolation == 'min':
                varValue[iCell] = np.min(var_adj)
            else:
                sys.exit("ERROR: wrong extrapolation scheme! Set option x as idw or min!")

            keepCellMaskNew[iCell] = 1

    print ("{0:8d} cells left for extrapolation in total {1:8d} cells".format(nCells-np.count_nonzero(keepCellMask),  nCells))
dataset.variables[var_name][0,:] = varValue # Put updated array back into file.
# === Clean-up =============

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
thiscommand = thiscommand+";  {} successfully extrapolated using the {} method".format(var_name, extrapolation)
if hasattr(dataset, 'history'):
   newhist = '\n'.join([thiscommand, getattr(dataset, 'history')])
else:
   newhist = thiscommand
setattr(dataset, 'history', newhist)


dataset.close()
