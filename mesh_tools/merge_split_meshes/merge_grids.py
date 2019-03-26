#!/usr/bin/env python
'''
Tool to merge 2 MPAS non-contiguous meshes together into a single file
'''

import sys
import numpy as np
import netCDF4
import argparse
import math
from collections import OrderedDict
import scipy.spatial
import time
from datetime import datetime


#print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n"
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.description = __doc__
parser.add_argument("-1", dest="file1", help="name of file 1", metavar="FILENAME")
parser.add_argument("-2", dest="file2", help="name of file 2", metavar="FILENAME")
parser.add_argument("-o", dest="outFile", help="name of output file", default="merged_mpas.nc", metavar="FILENAME")
#for option in parser.option_list:
#    if option.default != ("NO", "DEFAULT"):
#        option.help += (" " if option.help else "") + "[default: %default]"
options = parser.parse_args()



f1 = netCDF4.Dataset(options.file1)
nCells1 = len(f1.dimensions['nCells'])
nEdges1 = len(f1.dimensions['nEdges'])
nVertices1 = len(f1.dimensions['nVertices'])
Time1= len(f1.dimensions['Time'])

f2 = netCDF4.Dataset(options.file2)
nCells2 = len(f2.dimensions['nCells'])
nEdges2 = len(f2.dimensions['nEdges'])
nVertices2 = len(f2.dimensions['nVertices'])
Time2= len(f2.dimensions['Time'])

if Time1 != Time2:
  sys.exit("ERROR: The two files have different lengths of the Time dimension.")
if len(f1.dimensions['vertexDegree']) != len(f2.dimensions['vertexDegree']):
  sys.exit("ERROR: The two files have different lengths of the vertexDegree dimension.")
if len(f1.dimensions['nVertLevels']) != len(f2.dimensions['nVertLevels']):
  sys.exit("ERROR: The two files have different lengths of the nVertLevels dimension.")


# Create new file
fout = netCDF4.Dataset(options.outFile, "w", format="NETCDF3_CLASSIC")

# add merged dimensions
print("Adding merged dimensions to new file.")
fout.createDimension('nCells', nCells1+nCells2)
fout.createDimension('nEdges', nEdges1+nEdges2)
fout.createDimension('nVertices', nVertices1+nVertices2)
fout.createDimension('TWO', 2)
fout.createDimension('vertexDegree', len(f1.dimensions['vertexDegree']))
if 'StrLen' in f1.dimensions:
   fout.createDimension('StrLen', len(f1.dimensions['StrLen']))
maxEdges = max(len(f1.dimensions['maxEdges']), len(f2.dimensions['maxEdges']))
fout.createDimension('maxEdges', maxEdges)
fout.createDimension('maxEdges2', maxEdges*2)
fout.createDimension('nVertLevels', len(f1.dimensions['nVertLevels']))
fout.createDimension('nVertInterfaces', len(f1.dimensions['nVertInterfaces']))

fout.createDimension('Time', size=None) # make unlimited dimension


# compare list of variables
vars1 = f1.variables
vars2 = f2.variables

# only copy variables common to both files
for varname in vars1:
   if varname in vars2:
      print("Merging variable {}".format(varname))
      if f1.variables[varname].dimensions != f2.variables[varname].dimensions:
         sys.exit("Error: Variable {} has different dimensions in the two files.").format(varname)

      theVar = f1.variables[varname]
      newVar = fout.createVariable(varname, theVar.dtype, theVar.dimensions)
      # (Assuming here that nCells, nEdges, and nVertices are never both in a variable)
      # now assign value
      if 'nCells' in theVar.dimensions:
          ind = theVar.dimensions.index('nCells')
          tup1 = ()
          tup2 = ()
          tupMerge = ()
          for ind in range(len(theVar.dimensions)):
              if theVar.dimensions[ind] == 'nCells':
                  tup1 += (slice(0,nCells1),)
                  tup2 += (slice(0,nCells2),)
                  tupMerge += (slice(nCells1, nCells1+nCells2),)
              else:
                  tup1 += (slice(None),)
                  tup2 += (slice(None),)
                  tupMerge += (slice(None),)
          newVar[tup1] = f1.variables[varname][tup1]
          newVar[tupMerge] = f2.variables[varname][tup2]
      elif 'nEdges' in theVar.dimensions:
          ind = theVar.dimensions.index('nEdges')
          tup1 = ()
          tup2 = ()
          tupMerge = ()
          for ind in range(len(theVar.dimensions)):
              if theVar.dimensions[ind] == 'nEdges':
                  tup1 += (slice(0,nEdges1),)
                  tup2 += (slice(0,nEdges2),)
                  tupMerge += (slice(nEdges1, nEdges1+nEdges2),)
              else:
                  tup1 += (slice(None),)
                  tup2 += (slice(None),)
                  tupMerge += (slice(None),)
          newVar[tup1] = f1.variables[varname][tup1]
          newVar[tupMerge] = f2.variables[varname][tup2]
      elif 'nVertices' in theVar.dimensions:
          ind = theVar.dimensions.index('nVertices')
          tup1 = ()
          tup2 = ()
          tupMerge = ()
          for ind in range(len(theVar.dimensions)):
              if theVar.dimensions[ind] == 'nVertices':
                  tup1 += (slice(0,nVertices1),)
                  tup2 += (slice(0,nVertices2),)
                  tupMerge += (slice(nVertices1, nVertices1+nVertices2),)
              else:
                  tup1 += (slice(None),)
                  tup2 += (slice(None),)
                  tupMerge += (slice(None),)
          newVar[tup1] = f1.variables[varname][tup1]
          newVar[tupMerge] = f2.variables[varname][tup2]
      else:
          # just take file 1's version
          newVar[:] = theVar[:]

      # Indexes need adjusting:
      if varname == "indexToCellID":
          newVar[nCells1:] += nCells1
      elif varname == "indexToEdgeID":
          newVar[nEdges1:] += nEdges1
      elif varname == "indexToVertexID":
          newVar[nVertices1:] += nVertices1
      elif varname == "cellsOnEdge":
          part2 = newVar[nEdges1:,:]
          part2[part2>0] += nCells1
          newVar[nEdges1:,:] = part2
      elif varname == "edgesOnCell":
          part2 = newVar[nCells1:,:]
          part2[part2>0] += nEdges1
          newVar[nCells1:,:] = part2
      elif varname == "edgesOnEdge":
          part2 = newVar[nEdges1:,:]
          part2[part2>0] += nEdges1
          newVar[nEdges1:,:] = part2
      elif varname == "cellsOnCell":
          part2 = newVar[nCells1:,:]
          part2[part2>0] += nCells1
          newVar[nCells1:,:] = part2
      elif varname == "verticesOnCell":
          part2 = newVar[nCells1:,:]
          part2[part2>0] += nVertices1
          newVar[nCells1:,:] = part2
      elif varname == "verticesOnEdge":
          part2 = newVar[nEdges1:,:]
          part2[part2>0] += nVertices1
          newVar[nEdges1:,:] = part2
      elif varname == "edgesOnVertex":
          part2 = newVar[nVertices1:,:]
          part2[part2>0] += nEdges1
          newVar[nVertices1:,:] = part2
      elif varname == "cellsOnVertex":
          part2 = newVar[nVertices1:,:]
          part2[part2>0] += nCells1
          newVar[nVertices1:,:] = part2


# add some needed attributes
fout.on_a_sphere = "NO"
fout.sphere_radius = 0.0
fout.is_periodic = "NO"
# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
setattr(fout, 'history', thiscommand )
fout.close()
f1.close()
f2.close()

print('\nMerge completed.')

