#!/usr/bin/env python
"""
Script to extrapolate arbitrary MALI variable

Created on Mon Feb1 2021

@author: Matt Hoffman, Trevor Hillebrand
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import sys
from datetime import datetime


parser = OptionParser(
    description='Extrapolate a variable into undefined regions')

parser.add_option("-f", "--file", dest="nc_file",
                  help="the mpas file to write to")
parser.add_option("-v", "--variable", dest="var_name",
                  help="the MALI variable you want to extrapolate")
parser.add_option("-m", "--method", dest="extrapolation", default='min',
                  help="idw, min, or value method of extrapolation")
parser.add_option("-s", "--set_value", dest="set_value", default=None,
                  help=("value to set variable to outside "
                        "keepCellMask when using -v value"))

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

dataset = Dataset(options.nc_file, 'r+')
dataset.set_auto_mask(False)
var_name = options.var_name
extrapolation = options.extrapolation
# Extrapolation
nCells = len(dataset.dimensions['nCells'])
if 'thickness' in dataset.variables.keys():
    thickness = dataset.variables['thickness'][0,:]
    bed = dataset.variables["bedTopography"][0,:]
    geometry = True
else:
    geometry = False

cellsOnCell = dataset.variables['cellsOnCell'][:]
nEdgesOnCell = dataset.variables['nEdgesOnCell'][:]
xCell = dataset.variables["yCell"][:]
yCell = dataset.variables["xCell"][:]

if dataset.variables[var_name].ndim == 2:
    has_vertical_dim = False
    n_vert = 1
    varValue = dataset.variables[var_name][0, :]
elif dataset.variables[var_name].ndim == 3:
    has_vertical_dim = True
    print(dataset.variables[var_name].dimensions[2])
    vert_dim_name = dataset.variables[var_name].dimensions[2]
    n_vert = len(dataset.dimensions[vert_dim_name])
    print(f"This variable has a vertical dimension "
          f"of {vert_dim_name} with size {n_vert}")
    print("")
else:
    print("Unexpected number of dimension in variable")

for v in range(n_vert):
    if has_vertical_dim == True:
        print(f"Processing vertical level number {v}")
        varValue = dataset.variables[var_name][0, :, v]
    else:
        varValue = dataset.variables[var_name][0, :]

    # Define region of good data to extrapolate from.
    # Different methods for different variables
    if var_name in ["effectivePressure", "beta", "muFriction"]:
        groundedMask = (thickness > (-1028.0 / 910.0 * bed))
        keepCellMask = np.copy(groundedMask) * np.isfinite(varValue)

        for iCell in range(nCells):
            for n in range(nEdgesOnCell[iCell]):
                jCell = cellsOnCell[iCell, n] - 1
                if (groundedMask[jCell] == 1):
                    keepCellMask[iCell] = 1
                    continue
        # ensure zero muFriction does not get extrapolated
        keepCellMask *= (varValue > 0)
        # Get rid of invalid values
        keepCellMask *= (varValue < 1.e12)
    elif var_name in ["floatingBasalMassBal"]:
        if geometry:
            floatingMask = (thickness <= (-1028.0 / 910.0 * bed))
            keepCellMask = floatingMask * (varValue != 0.0)
        else:
            keepCellMask = (varValue != 0.0)
    else:
        keepCellMask = (thickness > 0.0)

    # make a copy to edit that will be used later
    keepCellMaskNew = np.copy(keepCellMask)
    # make a copy to edit that can be edited without changing the original
    keepCellMaskOrig = np.copy(keepCellMask)

    # recursive extrapolation steps:
    # 1) find cell A with mask = 0
    # 2) find how many surrounding cells have nonzero mask, and
    #    their indices (this will give us the cells from upstream)
    # 3) use the values for nonzero mask upstream
    #    cells toextrapolate the value for cell A
    # 4) change the mask for A from 0 to 1
    # 5) Update mask
    # 6) go to step 1)

    print("\nStart {} extrapolation using {} method".format(var_name, extrapolation))
    if extrapolation == 'value':
        varValue[np.where(np.logical_not(keepCellMask))] = float(options.set_value)
    else:
        while np.count_nonzero(keepCellMask) != nCells:
            keepCellMask = np.copy(keepCellMaskNew)
            searchCells = np.where(keepCellMask==0)[0]
            varValueOld = np.copy(varValue)

            for iCell in searchCells:
                neighborcellID = cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1
                # Important: ignore the phantom "neighbors" that are
                # off the edge of the mesh (0 values in cellsOnCell)
                neighborcellID = neighborcellID[neighborcellID>=0]
                # active cell mask
                mask_for_idx = keepCellMask[neighborcellID]
                mask_nonzero_idx, = np.nonzero(mask_for_idx)
                # id for nonzero beta cells
                nonzero_id = neighborcellID[mask_nonzero_idx]
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
                        sys.exit("ERROR: wrong extrapolation scheme!"
                                 " Set option m as idw or min!")

                    keepCellMaskNew[iCell] = 1

            print(f"{nCells-np.count_nonzero(keepCellMask)} cells left "
                  f"for extrapolation in total {nCells} cells")

    # Write updated array to file
    if has_vertical_dim == True:
        dataset.variables[var_name][0,:,v] = varValue
    else:
        dataset.variables[var_name][0,:] = varValue

# === Clean-up =============

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") \
              + ": " + " ".join(sys.argv[:])
thiscommand = thiscommand + f";  {var_name} successfully extrapolated " + \
              f"using the {extrapolation} method"
if hasattr(dataset, 'history'):
   newhist = '\n'.join([thiscommand, getattr(dataset, 'history')])
else:
   newhist = thiscommand
setattr(dataset, 'history', newhist)

dataset.close()
