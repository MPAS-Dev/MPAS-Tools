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


def extrapolate_into_mask(varValue, keepCellMask, cellNeighbors,
                          neighborValid, xCell, yCell, method):
    """Flood-fill varValue into cells where keepCellMask == 0.

    Vectorized, ring-by-ring equivalent of the original per-cell neighbor
    loop: each empty cell adjacent to the filled region takes the 'min' or
    inverse-distance-weighted ('idw') value of its already-filled neighbors
    (values frozen from the previous ring). Only the active frontier is
    processed each ring, so the cost scales with the number of cells rather
    than cells times rings. Cells not connected to any valid data are left
    unchanged (the original loop would never terminate in that case).
    """
    Nc = np.where(neighborValid, cellNeighbors, 0)  # safe gather index
    varValue = varValue.copy()
    filled = keepCellMask.astype(bool).copy()
    if method == 'idw':
        with np.errstate(divide='ignore'):
            dist = np.sqrt((xCell[:, None] - xCell[Nc]) ** 2 +
                           (yCell[:, None] - yCell[Nc]) ** 2)
            invDist = np.where(neighborValid,
                               1.0 / np.where(dist == 0, np.inf, dist), 0.0)
    elif method != 'min':
        sys.exit("ERROR: wrong extrapolation scheme! Set option m as idw or min!")
    newly = np.where(filled)[0]
    while newly.size:
        cand = np.unique(cellNeighbors[newly][neighborValid[newly]])
        cand = cand[~filled[cand]]
        if cand.size == 0:
            break
        nbFilled = neighborValid[cand] & filled[Nc[cand]]
        neighVals = varValue[Nc[cand]]
        if method == 'idw':
            w = np.where(nbFilled, invDist[cand], 0.0)
            varValue[cand] = np.sum(w * neighVals, axis=1) / np.sum(w, axis=1)
        else:
            varValue[cand] = np.where(nbFilled, neighVals, np.inf).min(axis=1)
        filled[cand] = True
        newly = cand
    return varValue


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

rho_ice = 910.
rho_ocean = 1028.
large_number = 1.e12

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

dataset = Dataset(options.nc_file, 'r+')
dataset.set_auto_mask(False)
var_name = options.var_name
extrapolation = options.extrapolation
nCells = len(dataset.dimensions['nCells'])
if all([geom_var in dataset.variables.keys() for
       geom_var in ["thickness", "bedTopography"]]):
    thickness = dataset.variables['thickness'][0, :]
    bed = dataset.variables["bedTopography"][0, :]
    geometry = True
else:
    geometry = False

cellsOnCell = dataset.variables['cellsOnCell'][:]
nEdgesOnCell = dataset.variables['nEdgesOnCell'][:]
xCell = dataset.variables["yCell"][:]
yCell = dataset.variables["xCell"][:]

# Precompute 0-based neighbor indices once (phantom neighbors -> -1).
cellNeighbors = cellsOnCell - 1
neighborValid = cellNeighbors >= 0

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
        if not geometry:
            raise Exception(f"Extrapolating {var_name} requires thickness "
                            "and bedTopography, but one or both of these "
                            f"are missing from {options.nc_file}")
        groundedMask = (thickness > (-1. * rho_ocean / rho_ice * bed))
        keepCellMask = np.copy(groundedMask) * np.isfinite(varValue)

        # Keep any cell that has a grounded neighbor (vectorized dilation over
        # the real edges of each cell).
        withinEdges = (np.arange(cellsOnCell.shape[1])[None, :] <
                       nEdgesOnCell[:, None])
        groundedNeighbor = (withinEdges &
                            (groundedMask[cellsOnCell - 1] == 1)).any(axis=1)
        keepCellMask[groundedNeighbor] = 1
        # ensure zero muFriction does not get extrapolated
        keepCellMask *= (varValue > 0)
        # Get rid of invalid values
        keepCellMask *= (varValue < large_number)
    elif var_name in ["floatingBasalMassBal"]:
        if geometry:
            floatingMask = (thickness <= (-1. * rho_ocean / rho_ice * bed))
            keepCellMask = floatingMask * (varValue != 0.0)
        else:
            keepCellMask = (varValue != 0.0)
    else:
        keepCellMask = (thickness > 0.0)

    # make a copy to edit that will be used later
    keepCellMaskNew = np.copy(keepCellMask)
    # make a copy to edit that can be edited without changing the original
    keepCellMaskOrig = np.copy(keepCellMask)

    # Flood-fill extrapolation: repeatedly fill empty cells adjacent to the
    # filled region using the 'min' or 'idw' value of their filled neighbors.

    print("\nStart {} extrapolation using {} method".format(var_name, extrapolation))
    if extrapolation == 'value':
        varValue[np.where(np.logical_not(keepCellMask))] = float(options.set_value)
    else:
        varValue = extrapolate_into_mask(
            varValue, keepCellMask, cellNeighbors, neighborValid,
            xCell, yCell, extrapolation)

    # Write updated array to file
    if has_vertical_dim == True:
        dataset.variables[var_name][0, :, v] = varValue
    else:
        dataset.variables[var_name][0, :] = varValue

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
