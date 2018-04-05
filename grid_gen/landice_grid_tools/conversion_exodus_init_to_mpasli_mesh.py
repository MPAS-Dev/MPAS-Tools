#!/usr/bin/env python
"""
Created on Tue Feb 13 23:50:20 2018

@author: Tong Zhang, Matt Hoffman
"""

import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import scipy.spatial as spt

parser = OptionParser(description='Read the basal friction data in the exo file and put them back to MPAS mesh. WARNING: Change the SEACAS library dir to your own path! A simple usage example: conversion_exodus_init_to_mpasli_mesh.py -e ./antarctica.exo -o target.nc -a ./mpas_cellID.ascii -v beta -k grd')
parser.add_option("-e", "--exo", dest="exo_file", help="the exo input file")
parser.add_option("-a", "--ascii", dest="id_file", help="the ascii global id input file")
parser.add_option("-m", "--method", dest="conversion_method", default="id", help="two options: id or coord. The id method is recommended. The coord method may fail at points where x or y = 0 while x_exodus or y_exodus is not")
parser.add_option("-k", "--mask", dest="mask_scheme", help="two options: all or grd. The all method is to mask cells with ice thickness > 0 as 1. The grd method masks grounded cells as 1.")
parser.add_option("-o", "--out", dest="nc_file", help="the mpas input/output file")
parser.add_option("-v", "--variable", dest="var_name", help="the mpas variable you want to convert from an exodus file")
parser.add_option("-x", "--extra", dest="extrapolation", default="min", help="Two options: idw and min. idw is the Inverse Distance Weighting method, and min is the method that uses the minimum value of the surrounding cells.")
parser.add_option("-i", "--iter", dest="smooth_iter_num", default="3", help="Maximum number for the recursive smoothing. A larger number means a more uniform smoothing field and more running time.")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

import sys
sys.path.append('/home/tzhang/Apps/seacas/lib')
from exodus import exodus
mpas_exodus_var_dic = {"beta":"basal_friction", "thickness":"ice_thickness", "stiffnessFactor":"stiffening_factor"}
# A mapping between mpas and exodus file. Add one if you need to manipulate a different variable!
ice_density = 910.0
ocean_density = 1028.0
useKDTree = False
# useKDTree is to find the surrounding nearest cells. Thus it only works for uniform meshes. But it's still useful for testing
# Normally, always use the cellsOnCell method. If you do want to use this method, set useKDTree to True 

exo_var_name = mpas_exodus_var_dic[options.var_name]

dataset = Dataset(options.nc_file, 'r+')
x = dataset.variables['xCell'][:]
y = dataset.variables['yCell'][:]

if useKDTree:
    xy = np.array([x,y])
    tree = spt.cKDTree(xy.T)

exo = exodus(options.exo_file)

stride = np.array(exo.get_global_variable_values('stride'))
ordering = np.array(exo.get_global_variable_values('ordering'))
# if ordering = 1, exo data is in the column-wise manner, stride is the vertical layer number
# if ordering = 0, exo data is in the layer-wise manner, stride is the node number each layer

data_exo = np.array(exo.get_node_variable_values(exo_var_name,1))
# read data from the exo file

xyz_exo = exo.get_coords()
x_exo = np.array(xyz_exo[0]) * 1000
y_exo = np.array(xyz_exo[1]) * 1000
# change the unit of the exo coord data from km to m. Be careful if it changes in the future

if ordering == 1.0:
    print "column wise pattern"
    layer_num = int(stride)
    data_exo_layer = data_exo[::layer_num]
    x_exo_layer = x_exo[::layer_num]
    y_exo_layer = y_exo[::layer_num]
elif ordering == 0.0:
    print "layer wise pattern"
    node_num = int(stride)
    data_exo_layer = data_exo[0:node_num+1]
    x_exo_layer = x_exo[0:node_num+1]
    y_exo_layer = y_exo[0:node_num+1]
    layer_num = int(len(data_exo)/node_num)
else:
    print "The ordering is probably wrong"
# slice the exo data to get the MPAS data

node_num_layer = len(x_exo_layer)

#dataset.variables[options.var_name][0,:] = 1.0e-0
# set beta value to some uniform value before we put new data in it

if (options.conversion_method == 'coord'):
    print "use coordinate method"
    for i in range(node_num_layer):
        index_x, = np.where(abs(x[:]-x_exo_layer[i])/(abs(x[:])+1e-10)<1e-3)
        index_y, = np.where(abs(y[:]-y_exo_layer[i])/(abs(y[:])+1e-10)<1e-3)
        index_intersect = list(set(index_x) & set(index_y))
        index = index_intersect[0]
        if options.var_name == "beta":
            dataset.variables[options.var_name][0,index] = np.exp(data_exo_layer[i]) * 1000.0
        else:
            dataset.variables[options.var_name][0,index] = (data_exo_layer[i])

        # The beta unit of the mpas mesh is Pa yr/m, not SI. 
        # This method may fail at the point where x or y = 0, while x_exo or y_exo is not

elif (options.conversion_method == 'id'):
    print "use global id method. Need a global id file"
    usefullCellID = np.loadtxt(options.id_file,dtype='i')
    usefullCellID_array = usefullCellID[1::]
    # The first number in the file is the total number. skip it

    if options.var_name == "beta":
        dataset.variables[options.var_name][0,usefullCellID_array-1] = np.exp(data_exo_layer) * 1000.0
    else:
        dataset.variables[options.var_name][0,usefullCellID_array-1] = np.exp(data_exo_layer)
        # We need convert fortran indice to python indice

else:
    sys.exit("wrong conversion method! Set option m as id or coord!")

print "Successful in converting data from Exodus to MPAS!"

nCells = len(dataset.dimensions['nCells'])
thickness = dataset.variables['thickness'][0,:]
bedrock = dataset.variables['bedTopography'][0,:]
cellsOnCell = dataset.variables['cellsOnCell'][:]
nEdgesOnCell = dataset.variables['nEdgesOnCell'][:]

keepCellMask = np.zeros((nCells,), dtype=np.int8)

if options.mask_scheme == 'grd':
    keepCellMask[thickness*ice_density/ocean_density + bedrock > 0.0] = 1
# find the mask for grounded ice region
elif options.mask_scheme == 'all':
    keepCellMask[thickness > 0.0] = 1
else:
    sys.exit("wrong masking scheme! Set option k as all or grd!")

keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that will be used later
keepCellMaskOld = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original

# recursive extrapolation steps:
# 1) find cell A with mask = 0 
# 2) find how many surrounding cells have nonzero mask, and their indices (this will give us the cells from upstream)
# 3) use the values for nonzero mask upstream cells to extrapolate the value for cell A
# 4) change the mask for A from 0 to 1
# 5) Update mask
# 6) go to step 1)

print "\nStart extrapolation!"

while np.count_nonzero(keepCellMask) != nCells:
    
    keepCellMask = np.copy(keepCellMaskNew)
    searchCells = np.where(keepCellMask==0)[0]

    if useKDTree:
        for iCell in searchCells:
            x_tmp = x[iCell]
            y_tmp = y[iCell]
            neareastPointNum = nEdgesOnCell[iCell]+1
            dist, idx = tree.query([x_tmp,y_tmp],k=neareastPointNum) # KD tree take account of [x_tmp,y_tmp] itself
            mask_for_idx = keepCellMask[idx]
            mask_nonzero_idx, = np.nonzero(mask_for_idx)

            nonzero_id = idx[mask_nonzero_idx]
            nonzero_num = np.count_nonzero(mask_for_idx)

            if nonzero_num > 0:
                dataset.variables[options.var_name][0,iCell] = sum(dataset.variables[options.var_name][0,nonzero_id])/nonzero_num
                keepCellMask[iCell] = 1

        print ("{0:8d} cells left for extrapolation in total {1:8d} cells".format(nCells-np.count_nonzero(keepCellMask),  nCells))


    else:
        for iCell in searchCells:
            neighborCellId = cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1

            mask_for_idx = keepCellMask[neighborCellId] # active cell mask
            mask_nonzero_idx, = np.nonzero(mask_for_idx)

            nonzero_id = neighborCellId[mask_nonzero_idx] # id for nonzero beta cells
            nonzero_num = np.count_nonzero(mask_for_idx)


            assert len(nonzero_id) == nonzero_num

            if nonzero_num > 0:
                x_i = x[iCell]
                y_i = y[iCell]
                x_adj = x[nonzero_id]
                y_adj = y[nonzero_id]
                ds = np.sqrt((x_i-x_adj)**2+(y_i-y_adj)**2)
                assert np.count_nonzero(ds)==len(ds)
                var_adj = dataset.variables[options.var_name][0,nonzero_id]
                if options.extrapolation == 'idw':
                    var_interp = 1.0/sum(1.0/ds)*sum(1.0/ds*var_adj)
                    dataset.variables[options.var_name][0,iCell] = var_interp
                elif options.extrapolation == 'min':
                    var_adj_min = min(var_adj)
                    dataset.variables[options.var_name][0,iCell] = var_adj_min
                else:
                    sys.exit("wrong extrapolation scheme! Set option x as idw or min!")

                keepCellMaskNew[iCell] = 1


        print ("{0:8d} cells left for extrapolation in total {1:8d} cells".format(nCells-np.count_nonzero(keepCellMask),  nCells))


print "\nStart idw smoothing for extrapolated field!"

iter_num = 0
while iter_num < int(options.smooth_iter_num):

    searchCells = np.where(keepCellMaskOld==0)[0]

    if useKDTree:
        for iCell in searchCells:
            x_tmp = x[iCell]
            y_tmp = y[iCell]
            neareastPointNum = nEdgesOnCell[iCell]+1
            dist, idx = tree.query([x_tmp,y_tmp],k=neareastPointNum) # KD tree take account of [x_tmp,y_tmp] itself
            mask_for_idx = keepCellMask[idx]
            mask_nonzero_idx, = np.nonzero(mask_for_idx)

            nonzero_id = idx[mask_nonzero_idx]
            nonzero_num = np.count_nonzero(mask_for_idx)

            assert nonzero_num > 0

            dataset.variables[options.var_name][0,iCell] = sum(dataset.variables[options.var_name][0,nonzero_id])/nonzero_num

        print ("{0:3d} smoothing in total {1:3s} iters".format(iter_num,  options.smooth_iter_num))


    else:
        for iCell in searchCells:
            neighborCellId = cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1
            nonzero_idx, = np.nonzero(neighborCellId+1)
            cell_nonzero_id = neighborCellId[nonzero_idx] # neighbor cell id

            mask_for_idx = keepCellMask[cell_nonzero_id] # active cell mask
            mask_nonzero_idx, = np.nonzero(mask_for_idx)

            nonzero_id = cell_nonzero_id[mask_nonzero_idx] # id for nonzero beta cells
            nonzero_num = np.count_nonzero(mask_for_idx)

            assert len(nonzero_id) == nonzero_num
            assert nonzero_num > 0

            x_i = x[iCell]
            y_i = y[iCell]
            x_adj = x[nonzero_id]
            y_adj = y[nonzero_id]
            ds = np.sqrt((x_i-x_adj)**2+(y_i-y_adj)**2)
            assert np.count_nonzero(ds)==len(ds)
            var_adj = dataset.variables[options.var_name][0,nonzero_id]
            var_interp = 1.0/sum(1.0/ds)*sum(1.0/ds*var_adj)
            dataset.variables[options.var_name][0,iCell] = var_interp

        print ("{0:3d} smoothing in total {1:3s} iters".format(iter_num,  options.smooth_iter_num))

    iter_num = iter_num + 1

if iter_num == 0:
    print "\nNo smoothing! Iter number is 0!"

print "\nExtrapolation and smoothing finished!"

dataset.close()
