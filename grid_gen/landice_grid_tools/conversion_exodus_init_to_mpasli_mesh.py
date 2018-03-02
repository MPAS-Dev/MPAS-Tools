#!/usr/bin/env python
"""
Created on Tue Feb 13 23:50:20 2018

@author: Tong Zhang
"""

import sys
sys.path.append('/home/tzhang/Apps/seacas/lib')
from exodus import exodus
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import scipy.spatial as spt

parser = OptionParser(description='Read the basal friction data in the exo file and put them back to MPAS mesh. Change the SEACAS library dir to your own path!')
parser.add_option("-e", "--exo", dest="exo_file", help="the exo input file")
parser.add_option("-a", "--ascii", dest="id_file", help="the ascii global id input file")
parser.add_option("-m", "--method", dest="conversion_method", help="two options: id or coord. The id method is recommended. The coord method may fail at points where x or y = 0 while x_exodus or y_exodus is not")
parser.add_option("-o", "--out", dest="nc_file", help="the mpas output file")
parser.add_option("-v", "--variable", dest="var_name", help="the mpas variable you want to convert from an exodus file")
parser.add_option("-x", "--extra", dest="extrapolation", help="set 1 if you want do extrapolation for surrounding buffer region, 0 for turning it down. The current extrapolation method is a simple average of active ice cell values. A distance inversed method is still needed for variable resolution meshes.")
options, args = parser.parse_args()

mpas_exodus_var_dic = {"beta":"basal_friction", "thickness":"ice_thickness", "enhance:enhancement"}
# A mapping between mpas and exodus file. Add one if you need to manipulate a different variable!

exo_var_name = mpas_exodus_var_dic[options.var_name]

dataset = Dataset(options.nc_file, 'r+', format="NETCDF4")
x = dataset.variables['xCell'][:]
y = dataset.variables['yCell'][:]
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
    #thickness_exo_layer = thickness_exo[::layer_num]
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

#dataset.variables[options.var_name][0,:] = 1.0e-10
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
        elif options.var_name == "enhance":
            for l in range(layer_num):
                dataset.variables[options.var_name][0,index,l] = data_exo_layer[i]
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
    elif options.var_name == "enhance"
        for l in range(layer_num):
            dataset.variables[options.var_name][0,usefullCellID_array-1,l] = data_exo_layer
    else:
        dataset.variables[options.var_name][0,usefullCellID_array-1] = (data_exo_layer)
        # We need convert fortran indice to python indice

else:
    sys.exit("wrong conversion method! Use id or coord only!")


nCells = len(dataset.dimensions['nCells'])
thickness = dataset.variables['thickness'][0,:]
bedrock = dataset.variables['bedTopography'][0,:]
cellsOnCell = dataset.variables['cellsOnCell'][:]
nEdgesOnCell = dataset.variables['nEdgesOnCell'][:]

ice_density = 910.0
ocean_density = 1028.0

keepCellMask = np.zeros((nCells,), dtype=np.int8)
keepCellMask[:] = 0
keepCellMask[thickness*ice_density/ocean_density + bedrock > 0.0] = 1
# find the mask for grounded ice region
#keepCellMask[thickness > 0.0] = 1

keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original

# recursive extrapolation steps:
# 1) find cell A with mask = 0 (ice_thickness = 0)
# 2) find six adjacent cells around A
# 3) find how many surrounding cells have nonzero mask, and their indices
# 4) use beta for nonzero mask surrounding cells to extrapolate the beta for cell A
# 5) change the mask for A from 0 to 1
# 6) Update mask
# 7) go to step 1)

if options.extrapolation == "1":

    while np.count_nonzero(keepCellMask) != nCells:

        # The following section use KD Tree to find the surrounding nearest cells. Thus it only works for uniform meshes. But it's still useful for testing
        # Normally, always use the method in the else section. If you do want to use this method, change "if 0" below to "if 1"
        if 0:
            for iCell in range(nCells):
                if keepCellMask[iCell] == 0:
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
                        keepCellMaskNew[iCell] = 1;

            keepCellMask = np.copy(keepCellMaskNew)
            print "%d cells left for extrapolation!" % (nCells-np.count_nonzero(keepCellMask))

        else:
            for iCell in range(nCells):
                if keepCellMask[iCell] == 0:
                    neighborCellId = cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1
                    nonzero_idx, = np.nonzero(neighborCellId+1)
                    cell_nonzero_id = neighborCellId[nonzero_idx] # neighbor cell id

                    mask_for_idx = keepCellMask[cell_nonzero_id] # active cell mask
                    mask_nonzero_idx, = np.nonzero(mask_for_idx)

                    nonzero_id = cell_nonzero_id[mask_nonzero_idx] # id for nonzero beta cells
                    nonzero_num = np.count_nonzero(mask_for_idx)

                    if len(nonzero_id) == nonzero_num:

                        if nonzero_num > 0:
                            dataset.variables[options.var_name][0,iCell] = sum(dataset.variables[options.var_name][0,nonzero_id])/nonzero_num
                            keepCellMaskNew[iCell] = 1;

            keepCellMask = np.copy(keepCellMaskNew)
            print "%d cells left for extrapolation!" % (nCells-np.count_nonzero(keepCellMask))


dataset.close()
    
print "Successfull in converting data from Exodus to MPAS!"
    
