#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 23:50:20 2018

@author: tong
"""

import sys
sys.path.append('/home/tzhang/Apps/seacas/lib')
#sys.path.append('/home/tzhang/Apps/exomerge')
from exodus import exodus
#import exomerge
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import scipy.spatial as spt

parser = OptionParser(epilog='read the basal friction data in the exo file and put them back to MPAS mesh')
parser.add_option("-e", "--exo", dest="exo_file", help="the exo input file")
parser.add_option("-a", "--ascii", dest="id_file", help="the ascii global id input file")
parser.add_option("-m", "--method", dest="conversion_method", help="two options: id or coord")
parser.add_option("-o", "--out", dest="nc_file", help="the mpas output file")
options, args = parser.parse_args()


dataset = Dataset(options.nc_file, 'r+', format="NETCDF4")
#beta = dataset.variables['beta'][0,:]
#thickness = dataset.variables['thickness'][0,:]
# read the old basal friction data in the input MPAS file
#allCellID = dataset.variables['indexToCellID'][:]
# read the cell ID numbers in the input MPAS file
x = dataset.variables['xCell'][:]
y = dataset.variables['yCell'][:]
xy = np.array([x,y])
tree = spt.cKDTree(xy.T)

#usefullCellID = np.loadtxt(options.id_file)
# read the targeted global id numbers from the ascii input file
#usefullCellID = [int(i) for i in usefullCellID[1::]]
# convert the number from float to int

#model = exomerge.import_model(options.exo_file)
exo = exodus(options.exo_file)
#model.output_global_variables('global_variables.csv')

#stride = model.global_variables['stride'][0]
stride = np.array(exo.get_global_variable_values('stride'))
#ordering = model.global_variables['ordering'][0]
ordering = np.array(exo.get_global_variable_values('ordering'))
# if ordering = 1, exo data is in the column-wise manner, stride is the vertical layer number
# if ordering = 0, exo data is in the layer-wise manner, stride is the node number each layer

#basal_friction = model.get_node_field_values("basal_friction")
beta_exo = np.array(exo.get_node_variable_values('basal_friction',1))
thickness_exo = np.array(exo.get_node_variable_values('ice_thickness',1))
# read the basal_friction data from the exo file
#ice_thickness = model.get_node_field_values("ice_thickness")
#ice_thickness = exo.get_node_variable_values('ice_thickness',1)[:]

xyz_exo = exo.get_coords()
x_exo = np.array(xyz_exo[0]) * 1000
y_exo = np.array(xyz_exo[1]) * 1000
# change the unit of the exo coord data from km to m. Be careful if it changes in the future

if ordering == 1.0:
    print "column wise pattern"
    layer_num = int(stride)
    beta_exo_layer = beta_exo[::layer_num]
    thickness_exo_layer = thickness_exo[::layer_num]
    x_exo_layer = x_exo[::layer_num]
    y_exo_layer = y_exo[::layer_num]
#    ice_thickness_MPAS = ice_thickness[::layer_num]
elif ordering == 0.0:
    print "layer wise pattern"
    node_num = int(stride)
    beta_exo_layer = beta_exo[0:node_num+1]
    x_exo_layer = x_exo[0:node_num+1]
    y_exo_layer = y_exo[0:node_num+1]
else:
    print "The ordering is probably wrong"
# slice the exo data to get the MPAS data

node_num_layer = len(x_exo_layer)

#index_array = np.zeros(node_num_layer)
#thickness_array = np.zeros(node_num_layer)

dataset.variables['beta'][0,:] = 1.0e-10

if (options.conversion_method == 'coord'):
    print "use coordinate method"
    for i in range(node_num_layer):
        index_x, = np.where(abs(x[:]-x_exo_layer[i])<1e2)
        index_y, = np.where(abs(y[:]-y_exo_layer[i])<1e2)
        # perhaps there is a better criteria for identifying coords?
        index_intersect = list(set(index_x) & set(index_y))
        index = index_intersect[0]
        #index_array[i] = index
        #thickness_array[i] = thickness_exo_layer[i] 

        #dataset.variables['beta'][0,index] = thickness_exo_layer[i]
        dataset.variables['beta'][0,index] = np.exp(beta_exo_layer[i])
        # The beta unit of the mpas mesh is kPa yr/m ?
elif (options.conversion_method == 'id'):
    print "use global id method. A global id file have to supplied by prescribing the id_file"
    usefullCellID = np.loadtxt(options.id_file,dtype='i')
    usefullCellID_array = usefullCellID[1::]
    # The first number in the file is the total number. skip it
    #dataset.variables['beta'][0,usefullCellID_array-1] = np.exp(beta_exo_layer)
    dataset.variables['beta'][0,usefullCellID_array-1] = np.exp(beta_exo_layer)
    # We need convert fortran indice to python indice

else:
    print "wrong conversion method! Use id or coord only!"


nCells = len(dataset.dimensions['nCells'])
cullCell = np.zeros((nCells, ), dtype=np.int8)  # Initialize to cull no cells
thickness = dataset.variables['thickness'][0,:]
#cellsOnCell = dataset.variables['cellsOnCell'][:]
#nEdgesOnCell = dataset.variables['nEdgesOnCell'][:]

keepCellMask = np.copy(cullCell[:])
keepCellMask[:] = 0
keepCellMask[thickness > 0.0] = 1


keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original

#for iCell in range(len(nEdgesOnCell)):
#    if keepCellMask[iCell] == 0:
#        for neighbor in cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1:  # the -1 converts from the fortran indexing in the variable to python indexing
#            if neighbor >= 0 and keepCellMask[neighbor] == 1:  # if any neighbors are already being kept on the old mask then keep this cell too.  This will get a lot of the ice cells, but they have already been assigned so they don't matter.  What we care about is getting cells that are not ice that have one or more ice neighbors.
#                keepCellMaskNew[iCell] = 1
#
#
#dataset.variables['beta'][0,:]=keepCellMaskNew


#cellsOnCell_new = np.copy(cellsOnCell)

# recursive extrapolation steps:
# 0) build DK tree using cell x y coordinates
# 1) find cell A with mask = 0 (ice_thickness = 0)
# 2) find six adjacent cells around A
# 3) find how many surrounding cells have nonzero mask, and their indices
# 4) use beta for nonzero mask surrounding cells to extrapolate the beta for cell A
# 5) change the mask for A from 0 to 1
# 6) Update mask
# 7) go to step 1)

while np.count_nonzero(keepCellMask) != nCells:


    for iCell in range(nCells):
        if keepCellMask[iCell] == 0:
            x_tmp = x[iCell]
            y_tmp = y[iCell]
            dist, idx = tree.query([x_tmp,y_tmp],k=6)
            mask_for_idx = keepCellMask[idx]
            mask_nonzero_idx, = np.nonzero(mask_for_idx)

            nonzero_id = idx[mask_nonzero_idx]
            nonzero_num = np.count_nonzero(mask_for_idx)

            if nonzero_num > 0:
                dataset.variables['beta'][0,iCell] = sum(dataset.variables['beta'][0,nonzero_id])/nonzero_num
                keepCellMaskNew[iCell] = 1;

    keepCellMask = np.copy(keepCellMaskNew)
    print "%d cells left for extrapolation!" % (nCells-np.count_nonzero(keepCellMask))




dataset.close()
    
print "Successfull in converting data from Exodus to MPAS!"
    
