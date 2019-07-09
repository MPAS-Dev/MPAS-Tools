#!/usr/bin/env python
'''
Script to remove repeated time entries in globalStats or other output file that occur
due to inexact restarts used in conjunction with the adpative timestepper.
Requires 'daysSinceStart' field is available. (Could be modified to use xtime instead)
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
from netCDF4 import Dataset
import numpy as np
import argparse


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--file", dest="file", help="File to be cleaned.", metavar="FILE", default="globalStats.nc")
args = parser.parse_args()

f = Dataset(args.file, 'r')
days = f.variables['daysSinceStart'][:]
nt = len(f.dimensions['Time'])

keepInd = np.zeros((nt,))
keepInd[0] = 1
prevMaxDay = days[0]
for i in range(1,nt):
  if days[i] > prevMaxDay:
     keepInd[i] = 1
     prevMaxDay = days[i]
  else:
     keepInd[i] = 0
print("Keeping {} indices out of {}".format(int(keepInd.sum()), nt))
keepList = np.nonzero(keepInd)[0]

if int(keepInd.sum())==nt:
   print("No cleaning required.")
   sys.exit()

# Copy all fields to a new file
fnameCleaned=args.file+".cleaned"
fileout = Dataset(fnameCleaned, 'w')

for name in f.ncattrs():
    setattr(fileout, name, getattr(f, name) )
    print('Copied global attribute  {} = {}'.format(name, getattr(f, name)))

if hasattr(fileout, 'history'):
   setattr(fileout, 'history', sys.argv[:] )

fileout.sync()

print("---- Copying dimensions from input file to output file ----")
for dim in f.dimensions.keys():
   print(dim)
   if dim == 'Time':
      dimvalue = None  # netCDF4 won't properly get this with the command below (you need to use the isunlimited method)
   else:
      dimvalue = len(f.dimensions[dim])
   fileout.createDimension(dim, dimvalue)
fileout.sync()

print("---- Copying variables from input file to output file ----")
for varname in f.variables:
   print(varname)
   thevar = f.variables[varname]
   newVar = fileout.createVariable(varname, thevar.dtype, thevar.dimensions)
   if 'Time' in f.variables[varname].dimensions:
      if 'Time' == f.variables[varname].dimensions[0]:
          if len(f.variables[varname].dimensions) == 1:
             newVar[:] = thevar[keepList]
          else:
             newVar[:] = thevar[keepList,:]
      else:
          sys.exit("Error: 'Time' is in dimension list for variable {}, but it is not the first dimension.  Script needs improving to handle this case.".format(varname))
   else:
      newVar[:] = thevar[:]
print("----")

fileout.close()
f.close()

print("Complete.  Cleaned output written to {}".format(fnameCleaned))
