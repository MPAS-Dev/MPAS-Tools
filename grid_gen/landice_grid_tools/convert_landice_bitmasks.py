#!/usr/bin/env python
'''
Script to convert landice bit mask into individual masks for each bit and save them to the netcdf file.
Converts any of cellMask, edgeMask, vertexMask present in file.
'''
import numpy
from netCDF4 import Dataset

from optparse import OptionParser

# Copied from src/core_landice/shared/mpas_li_mask.F
li_mask_ValueIce                   =  32  # Giving this the highest current value so it is obvious during visualization
li_mask_ValueDynamicIce            =   2
li_mask_ValueFloating              =   4
li_mask_ValueMargin                =   8  # This is the last cell with ice.
li_mask_ValueDynamicMargin         =  16  # This is the last dynamically active cell with ice
li_mask_ValueInitialIceExtent      =   1
li_mask_ValueAlbanyActive          =  64  # These are locations that Albany includes in its solution
li_mask_ValueAlbanyMarginNeighbor  = 128  # This the first cell beyond the last active albany cell
li_mask_ValueGroundingLine         = 256


masks = {'ice': li_mask_ValueIce,
         'dynamicIce': li_mask_ValueDynamicIce,
         'floating': li_mask_ValueFloating,
         'margin': li_mask_ValueMargin,
         'dynamicMargin': li_mask_ValueDynamicMargin,
         'initialExtent': li_mask_ValueInitialIceExtent,
         'albanyActive': li_mask_ValueAlbanyActive,
         'albanyMarginNeighbor': li_mask_ValueAlbanyMarginNeighbor,
         'groundingLine': li_mask_ValueGroundingLine
        }


print "** Gathering information."
parser = OptionParser()
parser.add_option("-f", "--filename", dest="filename", help="file to visualize; default: output.nc", default="output.nc", metavar="FILE")
options, args = parser.parse_args()


inFile = Dataset(options.filename, 'r+')
nTime = len(inFile.dimensions['Time'])


if 'cellMask' in inFile.variables:
   for maskName in masks:
     varName = "cellMask_"+maskName
     if varName in inFile.variables:
        newMaskVar = inFile.variables[varName]
     else:
        newMaskVar = inFile.createVariable(varName, 'i', ('Time','nCells'))
     for t in range(nTime):
        newMaskVar[t,:] = (inFile.variables['cellMask'][t,:] & masks[maskName]) / masks[maskName]
   inFile.sync()
   print "cellMask converted to individual masks."

if 'edgeMask' in inFile.variables:
   for maskName in masks:
     varName = "edgeMask_"+maskName
     if varName in inFile.variables:
        newMaskVar = inFile.variables[varName]
     else:
        newMaskVar = inFile.createVariable(varName, 'i', ('Time','nEdges'))
     for t in range(nTime):
        newMaskVar[t,:] = (inFile.variables['edgeMask'][t,:] & masks[maskName]) / masks[maskName]
   inFile.sync()
   print "edgeMask converted to individual masks."

if 'vertexMask' in inFile.variables:
   for maskName in masks:
     varName = "vertexMask_"+maskName
     if varName in inFile.variables:
        newMaskVar = inFile.variables[varName]
     else:
        newMaskVar = inFile.createVariable(varName, 'i', ('Time','nVertices'))
     for t in range(nTime):
        newMaskVar[t,:] = (inFile.variables['vertexMask'][t,:] & masks[maskName]) / masks[maskName]
   inFile.sync()
   print "vertexMask converted to individual masks."


inFile.close()
