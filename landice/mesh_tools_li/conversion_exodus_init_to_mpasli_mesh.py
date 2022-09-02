#!/usr/bin/env python
"""
Script to convert Albany-Land Ice output file in Exodus format to an MPAS-Land Ice format mesh.

See below for variables that are supported.

The script will use an environment variable SEACAS to find the Exodus python module.
To install the exodus python module, follow instructions here:
    https://github.com/gsjaardema/seacas
Note that you can follow the shortened instructions under "Exodus", and it is not necessary
to install the entire SEACAS set of libraries.

Created on Tue Feb 13 23:50:20 2018

@author: Tong Zhang, Matt Hoffman, Trevor Hillebrand
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
from scipy.interpolate import interp1d
import os
import sys
from datetime import datetime

parser = OptionParser(description='Convert data from exodus file to MPAS mesh. WARNING: Change the SEACAS library dir to your own path! A simple usage example: conversion_exodus_init_to_mpasli_mesh.py -e antarctica.exo -o target.nc -a ./mpas_cellID.ascii -v beta')
parser.add_option("-e", "--exo", dest="exo_file", help="the exo input file")
parser.add_option("-a", "--ascii", dest="id_file", help="the ascii global id input file")
parser.add_option("-o", "--out", dest="nc_file", help="the mpas file to write to")
parser.add_option("-v", "--variable", dest="var_name", help="the mpas variable(s) you want to convert from an exodus file. May be a single variable or multiple variables comma-separated (no spaces)")
parser.add_option("--dynamics", dest="dynamics", action="store_true", help="Use to convert ice dynamics fields: beta, muFriction, stiffnessFactor, uReconstructX/Y.  If stiffnessFactor was not included in the optimzation, it will be skipped.")
parser.add_option("--thermal", dest="thermal", action="store_true", help="Use to convert thermal fields: temperature, surfaceTemperature, basalTemperature.  Only use when the Albany optimization included the thermal solution.")
parser.add_option("--geometry", dest="geometry", action="store_true", help="Use to convert geometry fields: thickness and corresponding adjustment to bedTopography.  Only use when the Albany optimization included adjustment to the ice thickness.")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

# Get path to SEACAS module.  Two options:
#   environment variables SEACAS
#   or hard-coded default location
SEACAS_path = os.getenv('SEACAS')
if SEACAS_path == None:
   #sys.path.append('/home/tzhang/Apps/seacas/lib')
   #sys.path.append('/Users/trevorhillebrand/Documents/mpas/seacas/lib/')
   #sys.path.append('/Users/mhoffman/software/seacas/install/lib')
   sys.path.append('/global/project/projectdirs/piscees/nightlyCoriCDash/build/TrilinosInstall/lib')  # path on Cori
   # Note: There is currently not an operational installation on Badger or Grizzly
else:
   sys.path.append(SEACAS_path+'/lib')

try:
    from exodus import exodus
except ModuleNotFoundError:
    from exodus3 import exodus

# Map and copy Exodus data to MPAS data

# Create dictionary of variables that are supported by the script
mpas_exodus_var_dic = {"beta":"basal_friction_log",
                       "muFriction":"mu_log", \
                       "stiffnessFactor":"stiffening_factor", \
                       "uReconstructX":"solution_1", \
                       "uReconstructY":"solution_2", \
                       "temperature":"temperature", \
                       "basalTemperature":"temperature", \
                       "surfaceTemperature":"temperature", \
                       "thickness":"ice_thickness"}
                       #"surfaceAirTemperature":"surface_air_temperature", \ #use with caution!
# A mapping between mpas and exodus file. Add one if you need to manipulate a different variable!
ice_density = 910.0
ocean_density = 1028.0

if not options.nc_file:
   sys.exit("ERROR: an MPAS file for writing to was not specified with --out or -o.")
dataset = Dataset(options.nc_file, 'r+')
xCell = dataset.variables['xCell'][:]
yCell = dataset.variables['yCell'][:]
nCells = dataset.dimensions['nCells'].size

if not options.exo_file:
   sys.exit("ERROR: an Albany optimization exodus file was not specified with --exo or -e.")
exo = exodus(options.exo_file)

# Get Exodus coordinate arrays
xyz_exo = exo.get_coords()
x_exo = np.array(xyz_exo[0]) * 1000.0
y_exo = np.array(xyz_exo[1]) * 1000.0
# change the unit of the exo coord data from km to m. Be careful if it changes in the future

# Determine Exodus data ordering scheme
ordering = np.array(exo.get_global_variable_values('ordering'))

# Check that Albany and Exodus layer thicknesses are the same
exo_layer_name = []
exo_layer_thick = []
for glob_var_name in exo.get_global_variable_names():
    if "layer_thickness_ratio" in glob_var_name:
        exo_layer_name.append(glob_var_name)
        exo_layer_thick.append(exo.get_global_variable_values(glob_var_name))

# Albany layering is in reverse order from mpas
exo_layer_thick = np.flip(np.reshape(np.array(exo_layer_thick), [len(exo_layer_thick),]))
mpas_layer_thick = np.ma.filled(dataset.variables['layerThicknessFractions'][:])
if mpas_layer_thick.any() != exo_layer_thick.any():
    sys.exit("ERROR: Albany layer_thickness_ratio does not match MPAS layerThicknessFractions! Aborting")

# Read cellID_array from ascii file
if not options.id_file:
   sys.exit("ERROR: Cell ID file was not specified with --ascii or -a")
print("Reading global id file {}".format(options.id_file))
cellID = np.loadtxt(options.id_file,dtype='i')
# The first number in the file is the total number, which is the same as the stride
cellID_array = cellID[1::]
stride = cellID[0]  # Get stride from cell id file instead of from exo file

# Parse variable names from options
var_names = []
if options.dynamics:
    var_names.append("beta")
    var_names.append("muFriction")
    var_names.append("stiffnessFactor")
    var_names.append("uReconstructX")
    var_names.append("uReconstructY")
if options.thermal:
    var_names.append("temperature")
    var_names.append("surfaceTemperature")
    var_names.append("basalTemperature")
if options.geometry:
    var_names.append("thickness")

if len(var_names)>0 and options.var_name:
    sys.exit("ERROR: Specifying one or more of --dynamics --thermal --geometry and also specifying variables through --variable is not supported.")
if options.var_name:
    var_names = options.var_name.split(',') #parse variable names into a list

if len(var_names) == 0:
   sys.exit("ERROR: No variables for conversion were specified.  Specify variables with --dynamics --thermal --geometry or --variable.")

print("\nVariables to convert: ", *var_names)
vars_converted = []
vars_not_converted = []

# Loop through the variables, convert from Albany to MPAS, interpolate temperature
# from Albany to MPAS vertical layers, extrapolate and smooth beta and stiffnessFactor
for var_name in var_names:
    #set appropriate methods for smoothing and extrapolation
    if var_name == "beta":
        smooth_iter_num = 0
        extrapolation = "min"
    elif var_name == "muFriction":
        smooth_iter_num = 0
        extrapolation = "min"
    elif var_name == "stiffnessFactor":
        smooth_iter_num = 0
        extrapolation = "idw"
    else:
        smooth_iter_num = 0
        extrapolation = None

    print("\n---------------")
    if not var_name in mpas_exodus_var_dic:
        sys.exit("ERROR: variable '{}' not supported by this script.  Supported variables are: {}".format(var_name, mpas_exodus_var_dic.keys()))
    if mpas_exodus_var_dic[var_name] in exo.get_node_variable_names() and var_name in dataset.variables.keys():
        exo_var_name = mpas_exodus_var_dic[var_name]
        data_exo = np.array(exo.get_node_variable_values(exo_var_name,1))

        # Get number of vertical layers from mpas output file.
        if len(np.shape(dataset.variables[var_name])) == 3:
            if var_name == "temperature":
                nVert_max = np.shape(dataset.variables[var_name])[2] + 1 #albany temperature is on diferent vertical grid than MPAS, and has one extra layer
            else:
                nVert_max = np.shape(dataset.variables[var_name])[2]
        else:
            nVert_max = 1

        albanyTemperature = np.zeros((nCells, nVert_max)) #initialize albanyTemperature array to fill in below
        dataset.variables[var_name][:] = 0 # Fill variable with zeros in order to ensure proper values in void

        print("Begin {} conversion".format(var_name))
        # loop through nVertLevels (or nVertInterfaces) and convert variables from Albany to MPAS units
        for nVert in np.arange(0, nVert_max):
            print('Converting layer/level {} of {}'.format(nVert+1, nVert_max))

            #Albany has inverted layer/level ordering relative to MPAS.
            if var_name == "surfaceTemperature":
                nVert_albany = nVert_max - 1
            elif var_name == "basalTemperature":
                nVert_albany = 0  # This not needed, because it will be 0 for basalTemperature, but adding this case to make that assumption explicit
            else:
                nVert_albany = nVert_max - nVert - 1

            #Extract data from exodus file
            # if ordering = 1, exo data is in the column-wise manner, stride is the vertical layer number
            # if ordering = 0, exo data is in the layer-wise manner, stride is the node number each layer
            if ordering == 1.0:
                layer_num = int(stride)
                data_exo_layer = data_exo[nVert_albany::layer_num]
                x_exo_layer = x_exo[nVert_albany::layer_num]
                y_exo_layer = y_exo[nVert_albany::layer_num]
            elif ordering == 0.0:
                node_num = int(stride)
                data_exo_layer = data_exo[nVert_albany*node_num:(nVert_albany+1)*node_num]
                x_exo_layer = x_exo[nVert_albany*node_num:(nVert_albany+1)*node_num]
                y_exo_layer = y_exo[nVert_albany*node_num:(nVert_albany+1)*node_num]
                if len(data_exo) % node_num != 0:
                   sys.exit("ERROR: The total number of Exodus nodes cannot be divided evenly by the number of nodes in a level.")
                layer_num = len(data_exo)//node_num
            else:
                sys.exit("Invalid ordering in Exodus file.  Ordering must be 0 or 1.")

            node_num_layer = len(x_exo_layer)

            if var_name == "beta":
                dataset.variables[var_name][0,cellID_array-1] = np.exp(data_exo_layer) * 1000.0
            elif var_name == "muFriction":
                dataset.variables[var_name][0,cellID_array-1] = np.exp(data_exo_layer)
            elif var_name == "stiffnessFactor":
                dataset.variables[var_name][0,cellID_array-1] = np.exp(data_exo_layer)
            elif var_name == "uReconstructX" or var_name == "uReconstructY":
                dataset.variables[var_name][0,cellID_array-1, nVert] = data_exo_layer / (60. * 60. * 24 * 365)
            elif var_name == "thickness":
                print("WARNING: thickness conversion is still experimental!  Carefully check results before using.")
                # We have to be careful to not change MPAS geometry in the extended layer of cells - only touch it where the MPAS file already had nonzero thickness
                thkMask = dataset.variables[var_name][0,cellID_array-1] > 1.0
                dataset.variables[var_name][0,cellID_array-1] = data_exo_layer * 1000.0 * thkMask
                # change bedTopography also when we change thickness, if that field exists
                if 'bedTopography' in dataset.variables.keys():
                    thicknessOrig = np.copy(dataset.variables[var_name][0,cellID_array-1])
                    bedTopographyOrig = np.copy(dataset.variables['bedTopography'][0,cellID_array-1])
                    surfaceTopographyOrig = thicknessOrig + bedTopographyOrig
                    dataset.variables['bedTopography'][0,cellID_array-1] = surfaceTopographyOrig - data_exo_layer * 1000.0
            elif var_name == "temperature":
                albanyTemperature[cellID_array-1, nVert] = data_exo_layer #interpolate onto MPAS mesh in interpolateTemperature
            elif var_name == "surfaceTemperature":
                dataset.variables[var_name][0,cellID_array-1] = data_exo_layer
            elif var_name in dataset.variables.keys() == False:
                sys.exit("Unsupported variable requested for conversion.")
            else:
                dataset.variables[var_name][0,cellID_array-1] = data_exo_layer

        if var_name == "temperature":
            # Make generic equally-spaced layers for albany and MPAS. This works
            # even for unequally-spaced layers because this is a linear interpolation
            # and MPAS temperature layers are always halfway between Albany temperature
            # layers.
            albany_layers = np.arange(0, nVert_max)
            MPAS_layers = np.arange(1, nVert_max) - 0.5 # MPAS and albany temperature layers are staggered
            temperatureInterpolant = interp1d(albany_layers, albanyTemperature, axis=1)
            dataset.variables["temperature"][0,:,:] = temperatureInterpolant(MPAS_layers)

            # Add reasonable (non-zero) temperatures outside of ice mask. Make this a linear
            # interpolation between min(273.15, surfaceAirTemperature) at the surface and 268.15K at the bed.
            surfaceAirTemperature = dataset.variables["surfaceAirTemperature"][0,:]
            surfaceAirTemperature[surfaceAirTemperature > 273.15] = 273.15

            for nLayer in np.arange(0, dataset.dimensions["nVertLevels"].size):
                tempInterpCells = dataset.variables["temperature"][0,:,nLayer] == 0.0
                dataset.variables["temperature"][0,tempInterpCells,nLayer] = (1 - np.sum(dataset.variables["layerThicknessFractions"][0:nLayer+1])) * \
                                surfaceAirTemperature[tempInterpCells] + np.sum(dataset.variables["layerThicknessFractions"][0:nLayer+1]) * 268.15 

            print('\nTemperature interpolation complete')

        # Extrapolate and smooth beta and stiffnessFactor fields
        if var_name in ["beta", "muFriction", "stiffnessFactor"]:
            # Extrapolation
            nCells = len(dataset.dimensions['nCells'])
            thickness = dataset.variables['thickness'][0,:]
            cellsOnCell = dataset.variables['cellsOnCell'][:]
            nEdgesOnCell = dataset.variables['nEdgesOnCell'][:]
            if 'bedTopography' in dataset.variables.keys():
                bedrock = dataset.variables['bedTopography'][0,:]
            else:
                bedrock = np.zeros(np.shape(thickness))
            varValue = dataset.variables[var_name][0,:] # get variable array value so we can work from memory and not disk

            keepCellMask = np.zeros((nCells,), dtype=np.int8)

            keepCellMask[varValue > 0.0] = 1 # Define region to keep as anywhere the optimization returned a value

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

            # Smoothing
            iter_num = 0
            varValue = dataset.variables[var_name][0,:] # get variable array value so we can work from memory and not disk
            print("\nStart {} smoothing".format(var_name))
            while iter_num < int(smooth_iter_num):
                print("\nSmoothing iteration {} of {}".format(iter_num+1,  smooth_iter_num))
                searchCells = np.where(keepCellMaskOrig==0)[0]

                for iCell in searchCells:
                    neighborcellID = cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1
                    neighborcellID = neighborcellID[neighborcellID>=0] # Important: ignore the phantom "neighbors" that are off the edge of the mesh (0 values in cellsOnCell)
                    nonzero_idx, = np.nonzero(neighborcellID+1)
                    cell_nonzero_id = neighborcellID[nonzero_idx] # neighbor cell id

                    mask_for_idx = keepCellMask[cell_nonzero_id] # active cell mask
                    mask_nonzero_idx, = np.nonzero(mask_for_idx)

                    nonzero_id = cell_nonzero_id[mask_nonzero_idx] # id for nonzero beta cells
                    nonzero_num = np.count_nonzero(mask_for_idx)

                    assert len(nonzero_id) == nonzero_num
                    assert nonzero_num > 0

                    x_i = xCell[iCell]
                    y_i = yCell[iCell]
                    x_adj = xCell[nonzero_id]
                    y_adj = yCell[nonzero_id]
                    ds = np.sqrt((x_i-x_adj)**2+(y_i-y_adj)**2)
                    assert np.count_nonzero(ds)==len(ds)
                    var_adj = varValue[nonzero_id]
                    if extrapolation == 'idw':
                        var_interp = 1.0/sum(1.0/ds)*sum(1.0/ds*var_adj)
                    elif extrapolation == 'min':
                        var_interp = min(var_adj)
                    else:
                        sys.exit("ERROR: Smoothing is only for beta and stiffness for now. Set option i to 0 to disable smoothing!")
                    varValue[iCell] = var_interp

                iter_num = iter_num + 1
            dataset.variables[var_name][0,:] = varValue # Put updated array back into file.

            if iter_num == 0:
                print("\nNo smoothing. Iter number is 0.")

            print("\n{} extrapolation and smoothing finished".format(var_name))
        vars_converted.append(var_name)

    else:
        print("WARNING: No equivalent of variable {} found".format(var_name))
        vars_not_converted.append(var_name)
    print("---------------\n")

print("\nVariables successfully converted: ", *vars_converted)
print("Variables not converted: ", *vars_not_converted)
print("\nConversion, extrapolation, and smoothing complete. Enjoy!")

# === Clean-up =============

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
thiscommand = thiscommand+";  Variables successfully converted from Exodus to MPAS: " + ",".join(vars_converted)
if hasattr(dataset, 'history'):
   newhist = '\n'.join([thiscommand, getattr(dataset, 'history')])
else:
   newhist = thiscommand
setattr(dataset, 'history', newhist)


dataset.close()
