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
parser.add_option("-o", "--out", dest="nc_file", help="the mpas input/output file")
parser.add_option("-v", "--variable", dest="var_name", help="the mpas variable(s) you want to convert from an exodus file. May be 'all', a single variable, or multiple variables comma-separated (no spaces)")
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
   sys.path.append('/usr/projects/climate/SHARED_CLIMATE/software/badger/trilinos/2018-12-19/gcc-6.4.0/openmpi-2.1.2/lib')  # path on LANL Badger/Grizzly
else:
   sys.path.append(SEACAS_path+'/lib')

from exodus import exodus

# Map and copy Exodus data to MPAS data

# Create dictionary of variables that are supported by the script
mpas_exodus_var_dic = {"beta":"basal_friction", "thickness":"ice_thickness",\
                       "stiffnessFactor":"stiffening_factor", \
                       "basalTemperature":"temperature", \
                       "surfaceTemperature":"temperature", \
                       "temperature":"temperature", \
                       #"surfaceAirTemperature":"surface_air_temperature", \ #use with caution!
                       "uReconstructX":"solution_1", \
                       "uReconstructY":"solution_2"}
# A mapping between mpas and exodus file. Add one if you need to manipulate a different variable!
ice_density = 910.0
ocean_density = 1028.0

dataset = Dataset(options.nc_file, 'r+')
xCell = dataset.variables['xCell'][:]
yCell = dataset.variables['yCell'][:]
nCells = dataset.dimensions['nCells'].size

exo = exodus(options.exo_file)

stride = np.array(exo.get_global_variable_values('stride'))

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
    sys.exit("Albany layer_thickness_ratio does not match MPAS layerThicknessFractions! Aborting")

# Read cellID_array from ascii file
print("Reading global id file {}".format(options.id_file))
cellID = np.loadtxt(options.id_file,dtype='i')
cellID_array = cellID[1::]
# The first number in the file is the total number. skip it

# Parse variable names from options
var_names = []
if options.var_name == "all":
    for mpas_name in mpas_exodus_var_dic:
        var_names.append(mpas_name)
else:
    var_names = options.var_name.split(',') #parse variable names into a list

# Loop through the variables, convert from Albany to MPAS, interpolate temperature
# from Albany to MPAS vertical layers, extrapolate and smooth beta and stiffnessFactor
for var_name in var_names:
    #set appropriate methods for smoothing and extrapolation
    if var_name == "beta":
        smooth_iter_num = 3
        mask_scheme = "grd"
        extrapolation = "min"
    elif var_name == "stiffnessFactor":
        smooth_iter_num = 3
        mask_scheme = "all"
        extrapolation = "idw"
    else:
        smooth_iter_num = 0
        mask_scheme = "all"
        extrapolation = None

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

        print("\nBegin {} conversion".format(var_name))
        # loop through nVertLevels (or nVertInterfaces) and convert variables from Albany to MPAS units
        for nVert in np.arange(0, nVert_max):
            print('Converting layer/level {} of {}'.format(nVert+1, nVert_max))
            
            #Albany has inverted layer/level ordering relative to MPAS.
            if var_name == "surfaceTemperature":
                nVert_albany = int(stride)-1
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
                data_exo_layer = data_exo[0:node_num+1]
                x_exo_layer = x_exo[0:node_num+1]
                y_exo_layer = y_exo[0:node_num+1]
                layer_num = len(data_exo)//node_num
            else:
                sys.exit("Invalid ordering in Exodus file.  Ordering must be 0 or 1.")

            node_num_layer = len(x_exo_layer)

            if var_name == "beta":
                dataset.variables[var_name][0,cellID_array-1] = np.exp(data_exo_layer) * 1000.0
            elif var_name == "uReconstructX" or var_name == "uReconstructY":
                dataset.variables[var_name][0,cellID_array-1, nVert] = data_exo_layer / (60. * 60. * 24 * 365)
            elif var_name == "thickness":
                # change bedTopography also when we change thickness, if that field exists
                dataset.variables[var_name][0,cellID_array-1] = data_exo_layer * 1000.0
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
            albany_layers = np.arange(0, nVert_max)
            MPAS_layers = np.arange(1, nVert_max) - 0.5 # MPAS and albany temperature layers are staggered
            temperatureInterpolant = interp1d(albany_layers, albanyTemperature, axis=1)
            dataset.variables["temperature"][0,:,:] = temperatureInterpolant(MPAS_layers)
            print('\nTemperature interpolation complete')
        
        # Extrapolate and smooth beta and stiffnessFactor fields
        if var_name in ["beta", "stiffnessFactor"]:
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

            if mask_scheme == 'grd':
                keepCellMask[(thickness*ice_density/ocean_density + bedrock > 0.0) * (thickness>0.0)] = 1
            # find the mask for grounded ice region
            elif mask_scheme == 'all':
                keepCellMask[thickness > 0.0] = 1

            keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that will be used later
            keepCellMaskOrig = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original


            # recursive extrapolation steps:
            # 1) find cell A with mask = 0
            # 2) find how many surrounding cells have nonzero mask, and their indices (this will give us the cells from upstream)
            # 3) use the values for nonzero mask upstream cells to extrapolate the value for cell A
            # 4) change the mask for A from 0 to 1
            # 5) Update mask
            # 6) go to step 1)

            print("\nStart {} extrapolation using {} method!".format(var_name, extrapolation))
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
                            sys.exit("wrong extrapolation scheme! Set option x as idw or min!")

                        keepCellMaskNew[iCell] = 1

                print ("{0:8d} cells left for extrapolation in total {1:8d} cells".format(nCells-np.count_nonzero(keepCellMask),  nCells))
            dataset.variables[var_name][0,:] = varValue # Put updated array back into file.

            # Smoothing
            iter_num = 0
            varValue = dataset.variables[var_name][0,:] # get variable array value so we can work from memory and not disk
            print("\nStart {} smoothing!".format(var_name))
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
                        sys.exit("Smoothing is only for beta and stiffness for now. Set option i to 0 to disable smoothing!")
                    varValue[iCell] = var_interp

                iter_num = iter_num + 1
            dataset.variables[var_name][0,:] = varValue # Put updated array back into file.

            if iter_num == 0:
                print("\nNo smoothing! Iter number is 0!")

            print("\n{} extrapolation and smoothing finished!".format(var_name))

    else:
        print("\nNo equivalent of variable {} found".format(var_name))

print("\nConversion, extrapolation, and smoothing complete. Enjoy!")

# === Clean-up =============

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
if hasattr(dataset, 'history'):
   newhist = '\n'.join([thiscommand, getattr(dataset, 'history')])
else:
   newhist = thiscommand
setattr(dataset, 'history', newhist)


dataset.close()
