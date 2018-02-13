#!/usr/bin/env python

import sys
sys.path.append('/home/tzhang/Apps/seacas/lib')
sys.path.append('/home/tzhang/Apps/exomerge')
from exodus import exodus
import exomerge
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser

parser = OptionParser(epilog='read the basal friction data in the exo file and put them back to MPAS mesh')
parser.add_option("-e", "--exo", dest="exo_file", help="the exo input file")
parser.add_option("-a", "--ascii", dest="id_file", help="the ascii global id input file")
parser.add_option("-o", "--out", dest="nc_file", help="the mpas output file")
options, args = parser.parse_args()


dataset = Dataset(options.nc_file, 'r+', format="NETCDF4")
beta = dataset.variables['beta']
# read the old basal friction data in the input MPAS file
allCellID = dataset.variables['indexToCellID']
# read the cell ID numbers in the input MPAS file

usefullCellID = np.loadtxt(options.id_file)
# read the targeted global id numbers from the ascii input file
usefullCellID = [int(i) for i in usefullCellID[1::]]
# convert the number from float to int

model = exomerge.import_model(options.exo_file)
#model.output_global_variables('global_variables.csv')

stride = model.global_variables['stride'][0]
ordering = model.global_variables['ordering'][0]
# if ordering = 1, exo data is in the column-wise manner, stride is the vertical layer number
# if ordering = 0, exo data is in the layer-wise manner, stride is the node number each layer

basal_friction = model.get_node_field_values("basal_friction")
# read the basal_friction data from the exo file
#ice_thickness = model.get_node_field_values("ice_thickness")

if ordering == 1.0:
    layer_num = int(stride)
    basal_friction_MPAS = basal_friction[::layer_num]
elif ordering == 0.0:
    node_num = int(stride)
    basal_friction_MPAS = basal_friction[0:node_num+1]
else:
    print "The ordering is probably wrong"
# slice the exo data to get the MPAS data

basal_friction_MPAS_array = np.array(basal_friction_MPAS)
usefullCellID_array = np.array(usefullCellID)
# a conversion from list to array

dataset.variables['beta'][0,usefullCellID_array] = np.exp(basal_friction_MPAS_array) * (365.0*24.0*3600.0)
# modify the basal friction values at the targeted global ids, and change the unit to Pa s m^{-1}

dataset.close()

print "Successfull in converting data from Exodus to MPAS!"

