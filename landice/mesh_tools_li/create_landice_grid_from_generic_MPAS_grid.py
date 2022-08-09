#!/usr/bin/env python
"""
Script to create a grid with land ice variables from an MPAS grid.
Currently variable attributes are not copied.
This script could be modified to copy them (looping over dir(var), skipping over variable function names "assignValue", "getValue", "typecode").
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import sys, numpy
from netCDF4 import Dataset
from optparse import OptionParser
from datetime import datetime


sphere_radius = 6.37122e6 # earth radius, if needed

print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser()
parser.add_option("-i", "--in", dest="fileinName", help="input filename.  Defaults to 'grid.nc'", metavar="FILENAME")
parser.add_option("-o", "--out", dest="fileoutName", help="output filename.  Defaults to 'landice_grid.nc'", metavar="FILENAME")
parser.add_option("-l", "--level", dest="levels", help="Number of vertical levels to use in the output file.  Defaults to the number in the input file", metavar="FILENAME")
parser.add_option("-v", "--vert", dest="vertMethod", help="Method of vertical layer spacing: uniform, glimmer.  Glimmer spacing follows Eq. 35 of Rutt, I. C., M. Hagdorn, N. R. J. Hulton, and A. J. Payne (2009), The Glimmer community ice sheet model, J. Geophys. Res., 114, F02004, doi:10.1029/2008JF001015", default='glimmer', metavar="FILENAME")
parser.add_option("--beta", dest="beta", action="store_true", help="DEPRECATED")
parser.add_option("--effecpress", dest="effecpress", action="store_true", help="DEPRECATED")
parser.add_option("--diri", dest="dirichlet", action="store_true", help="Use this flag to include the fields 'dirichletVelocityMask', 'uReconstructX', 'uReconstructY' needed for specifying Dirichlet velocity boundary conditions in the resulting file.")
parser.add_option("--thermal", dest="thermal", action="store_true", help="Use this flag to include the fields 'temperature', 'surfaceAirTemperature', 'basalHeatFlux' needed for specifying thermal initial conditions in the resulting file.")
parser.add_option("--hydro", dest="hydro", action="store_true", help="Use this flag to include the fields 'waterThickness', 'tillWaterThickness', 'basalMeltInput', 'externalWaterInput', 'frictionAngle', 'waterPressure', 'waterFluxMask' needed for specifying hydro initial conditions in the resulting file.")
parser.add_option("--obs", dest="obs", action="store_true", help="Use this flag to include the observational fields observedSurfaceVelocityX, observedSurfaceVelocityY, observedSurfaceVelocityUncertainty, observedThicknessTendency, observedThicknessTendencyUncertainty, thicknessUncertainty needed doing optimizations constrained by obs velocities.")
options, args = parser.parse_args()

if not options.fileinName:
    print("No input filename specified, so using 'grid.nc'.")
    options.fileinName = 'grid.nc'
else:
    print("Input file is: {}".format(options.fileinName))
if not options.fileoutName:
    print("No output filename specified, so using 'landice_grid.nc'.")
    options.fileoutName = 'landice_grid.nc'
print('') # make a space in stdout before further output

# Get the input file
filein = Dataset(options.fileinName,'r')

# Define the new file to be output
fileout = Dataset(options.fileoutName,"w",format=filein.file_format)


# ============================================
# Copy over all of the netcdf global attributes
# ============================================
# Do this first as doing it last is slow for big files since adding
# attributes forces the contents to get reorganized.
print("---- Copying global attributes from input file to output file ----")
for name in filein.ncattrs():
  # sphere radius needs to be set to that of the earth if on a sphere
  if name == 'sphere_radius' and getattr(filein, 'on_a_sphere') == "YES             ":
    setattr(fileout, 'sphere_radius', sphere_radius)
    print('Set global attribute   sphere_radius = {}'.format(sphere_radius))
  elif name =='history':
    # Update history attribute of netCDF file
    newhist = '\n'.join([getattr(filein, 'history'), ' '.join(sys.argv[:]) ] )
    setattr(fileout, 'history', newhist )
  else:
    # Otherwise simply copy the attr
    setattr(fileout, name, getattr(filein, name) )
    print('Copied global attribute  {} = {}'.format(name, getattr(filein, name)))

# Update history attribute of netCDF file if we didn't above
if not hasattr(fileout, 'history'):
   setattr(fileout, 'history', sys.argv[:] )

fileout.sync()
print('') # make a space in stdout before further output


# ============================================
# Copy over all the dimensions to the new file
# ============================================
# Note: looping over dimensions seems to result in them being written in seemingly random order.
#       I don't think this matters but it is not aesthetically pleasing.
#       It may be better to list them explicitly as I do for the grid variables,
#       but this way ensures they all get included and is easier.
# Note: The UNLIMITED time dimension will return a dimension value of None with Scientific.IO.  This is what is supposed to happen.  See below for how to deal with assigning values to a variable with a unlimited dimension.  Special handling is needed with the netCDF module.
print("---- Copying dimensions from input file to output file ----")
for dim in filein.dimensions.keys():
    if dim == 'nTracers':
        pass  # Do nothing - we don't want this dimension
    elif (dim == 'nVertInterfaces'):
        pass  # Do nothing - this dimension will be handled below
    else:    # Copy over all other dimensions
      if dim == 'Time':
         dimvalue = None  # netCDF4 won't properly get this with the command below (you need to use the isunlimited method)
      elif (dim == 'nVertLevels'):
        if options.levels is None:
          # If nVertLevels is in the input file, and a value for it was not
          # specified on the command line, then use the value from the file (do nothing here)
          print("Using nVertLevels from the intput file: {}".format(len(filein.dimensions[dim])))
          dimvalue = len(filein.dimensions[dim])
        else:
          # if nVertLevels is in the input file, but a value WAS specified
          # on the command line, then use the command line value
          print("Using nVertLevels specified on the command line: {}".format(int(options.levels)))
          dimvalue = int(options.levels)
      else:
         dimvalue = len(filein.dimensions[dim])
      fileout.createDimension(dim, dimvalue)
# There may be input files that do not have nVertLevels specified, in which case
# it has not been added to the output file yet.  Treat those here.
if 'nVertLevels' not in fileout.dimensions:
   if options.levels is None:
       print("nVertLevels not in input file and not specified.  Using default value of 10.")
       fileout.createDimension('nVertLevels', 10)
   else:
       print("Using nVertLevels specified on the command line: {}".format(int(options.levels)))
       fileout.createDimension('nVertLevels', int(options.levels))
# Also create the nVertInterfaces dimension, even if none of the variables require it.
fileout.createDimension('nVertInterfaces', len(fileout.dimensions['nVertLevels']) + 1)  # nVertInterfaces = nVertLevels + 1
print('Added new dimension nVertInterfaces to output file with value of {}.'.format(len(fileout.dimensions['nVertInterfaces'])))

fileout.sync()
print('Finished creating dimensions in output file.\n') # include an extra blank line here

# ============================================
# Copy over all of the required grid variables to the new file
# ============================================
print("Beginning to copy mesh variables to output file.")
vars2copy = ['latCell', 'lonCell', 'xCell', 'yCell', 'zCell', 'indexToCellID', 'latEdge', 'lonEdge', 'xEdge', 'yEdge', 'zEdge', 'indexToEdgeID', 'latVertex', 'lonVertex', 'xVertex', 'yVertex', 'zVertex', 'indexToVertexID', 'cellsOnEdge', 'nEdgesOnCell', 'nEdgesOnEdge', 'edgesOnCell', 'edgesOnEdge', 'weightsOnEdge', 'dvEdge', 'dcEdge', 'angleEdge', 'areaCell', 'areaTriangle', 'cellsOnCell', 'verticesOnCell', 'verticesOnEdge', 'edgesOnVertex', 'cellsOnVertex', 'kiteAreasOnVertex']
# Add these optional fields if they exist in the input file
for optionalVar in ['meshDensity', 'gridSpacing', 'cellQuality', 'triangleQuality', 'triangleAngleQuality', 'obtuseTriangle']:
   if optionalVar in filein.variables:
      vars2copy.append(optionalVar)

for varname in vars2copy:
   print("- ", end='')
print("|")
for varname in vars2copy:
   thevar = filein.variables[varname]
   datatype = thevar.dtype
   newVar = fileout.createVariable(varname, datatype, thevar.dimensions)
   if filein.on_a_sphere == "YES             ":
     if varname in ('xCell', 'yCell', 'zCell', 'xEdge', 'yEdge', 'zEdge', 'xVertex', 'yVertex', 'zVertex', 'dvEdge', 'dcEdge'):
       newVar[:] = thevar[:] * sphere_radius / filein.sphere_radius
     elif varname in ('areaCell', 'areaTriangle', 'kiteAreasOnVertex'):
       newVar[:] = thevar[:] * (sphere_radius / filein.sphere_radius)**2
     else:
       newVar[:] = thevar[:]
   else: # not on a sphere
     newVar[:] = thevar[:]
   del newVar, thevar
   sys.stdout.write("* "); sys.stdout.flush()
fileout.sync()
print("|")
print("Finished copying mesh variables to output file.\n")

# ============================================
# Create the land ice variables (all the shallow water vars in the input file can be ignored)
# ============================================
nVertLevels = len(fileout.dimensions['nVertLevels'])
datatype = filein.variables['xCell'].dtype  # Get the datatype for double precision float
datatypeInt = filein.variables['indexToCellID'].dtype  # Get the datatype for integers
#  Note: it may be necessary to make sure the Time dimension has size 1, rather than the 0 it defaults to.  For now, letting it be 0 which seems to be fine.

# layerThicknessFractions
layerThicknessFractions = fileout.createVariable('layerThicknessFractions', datatype, ('nVertLevels', ))
layerThicknessFractionsData = numpy.zeros(layerThicknessFractions.shape)
# Assign default values to layerThicknessFractions.  By default they will be uniform fractions.  Users can modify them in a subsequent step, but doing this here ensures the most likely values are already assigned. (Useful for e.g. setting up Greenland where the state variables are copied over but the grid variables are not modified.)
if options.vertMethod == 'uniform':
   layerThicknessFractionsData[:] = 1.0 / nVertLevels
elif options.vertMethod == 'glimmer':
   nInterfaces = nVertLevels + 1
   layerInterfaces = numpy.zeros((nInterfaces,))
   for k in range(nInterfaces):
      layerInterfaces[k] = 4.0/3.0 * (1.0 - ((k+1.0-1.0)/(nInterfaces-1.0) + 1.0)**-2)
   for k in range(nVertLevels):
      layerThicknessFractionsData[k] = layerInterfaces[k+1] - layerInterfaces[k]
   print("Setting layerThicknessFractions to: {}".format(layerThicknessFractionsData))
else:
   sys.exit('Unknown method for vertical spacing method (--vert): '+options.vertMethod)

# explictly specify layer fractions
#layerThicknessFractionsData[:] = [0.1663,0.1516,0.1368,0.1221,0.1074,0.0926,0.0779,0.0632,0.0484,0.0337]

layerThicknessFractions[:] = layerThicknessFractionsData[:]


# With Scientific.IO.netCDF, entries are appended along the unlimited dimension one at a time by assigning to a slice.
# Therefore we need to assign to time level 0, and what we need to assign is a zeros array that is the shape of the new variable, exluding the time dimension!
newvar = fileout.createVariable('thickness', datatype, ('Time', 'nCells'))
newvar[0,:] = numpy.zeros( newvar.shape[1:] )
# These landice variables are stored in the mesh currently, and therefore do not have a time dimension.
#    It may make sense to eventually move them to state.
newvar = fileout.createVariable('bedTopography', datatype, ('Time', 'nCells'))
newvar[:] = numpy.zeros(newvar.shape)
newvar = fileout.createVariable('sfcMassBal', datatype, ('Time', 'nCells'))
newvar[:] = numpy.zeros(newvar.shape)
newvar = fileout.createVariable('floatingBasalMassBal', datatype, ('Time', 'nCells'))
newvar[:] = numpy.zeros(newvar.shape)
print('Added default variables: thickness, temperature, bedTopography, sfcMassBal, floatingBasalMassBal')

newvar = fileout.createVariable('beta', datatype, ('Time', 'nCells'))
newvar[:] = 1.0e8  # Give a default beta that won't have much sliding.
print('Added variable: beta')

newvar = fileout.createVariable('muFriction', datatype, ('Time', 'nCells'))
newvar[:] = 1.0e8  # Give a default mu that won't have much sliding.
print('Added variable: muFriction')

newvar = fileout.createVariable('effectivePressure', datatype, ('Time', 'nCells'))
newvar[:] = 1.0  # Give a default effective pressure of 1.0 so that, for the linear sliding law, beta = mu*effecpress = mu.
print('Added variable: effectivePressure')

newvar = fileout.createVariable('stiffnessFactor', datatype, ('Time', 'nCells'))
newvar[:] = 1.0  # Give default value
print('Added variable: stiffnessFactor')

newvar = fileout.createVariable('eigencalvingParameter', datatype, ('Time', 'nCells'))
newvar[:] = 3.14e16  # Give default value for eigencalvingParameter
print('Added variable: eigencalvingParameter')

newvar = fileout.createVariable('groundedVonMisesThresholdStress', datatype, ('Time', 'nCells'))
newvar[:] = 1.0e6  # Give default value
print('Added variable: groundedVonMisesThresholdStress')

newvar = fileout.createVariable('floatingVonMisesThresholdStress', datatype, ('Time', 'nCells'))
newvar[:] = 1.0e6  # Give default value
print('Added variable: floatingVonMisesThresholdStress')

newvar = fileout.createVariable('iceMask', datatype, ('Time', 'nCells'))
newvar[:] = 0
print('Added variable: iceMask')

if options.dirichlet:
   newvar = fileout.createVariable('dirichletVelocityMask', datatypeInt, ('Time', 'nCells', 'nVertInterfaces'))
   newvar[:] = 0  # default: no Dirichlet b.c.
   newvar = fileout.createVariable('uReconstructX', datatype, ('Time', 'nCells', 'nVertInterfaces',))
   newvar[:] = 0.0
   newvar = fileout.createVariable('uReconstructY', datatype, ('Time', 'nCells', 'nVertInterfaces',))
   newvar[:] = 0.0
   print('Added optional dirichlet variables: dirichletVelocityMask, uReconstructX, uReconstructY')

if options.thermal:
   newvar = fileout.createVariable('temperature', datatype, ('Time', 'nCells', 'nVertLevels'))
   newvar[:] = 273.15 # Give default value for temperate ice
   newvar = fileout.createVariable('surfaceAirTemperature', datatype, ('Time', 'nCells'))
   newvar[:] = 273.15 # Give default value for temperate ice
   newvar = fileout.createVariable('basalHeatFlux', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0 # Default to none (W/m2)
   print('Added optional thermal variables: temperature, surfaceAirTemperature, basalHeatFlux')

if options.hydro:
   newvar = fileout.createVariable('waterThickness', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('tillWaterThickness', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('basalMeltInput', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('externalWaterInput', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('frictionAngle', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('waterPressure', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('waterFluxMask', 'i', ('Time', 'nEdges'))
   newvar[:] = 0.0
   print('Added optional hydro variables: waterThickness, tillWaterThickness, meltInput, frictionAngle, waterPressure, waterFluxMask')

if options.obs:
   newvar = fileout.createVariable('observedSurfaceVelocityX', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('observedSurfaceVelocityY', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('observedSurfaceVelocityUncertainty', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('observedThicknessTendency', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('observedThicknessTendencyUncertainty', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   newvar = fileout.createVariable('thicknessUncertainty', datatype, ('Time', 'nCells'))
   newvar[:] = 0.0
   print('Added optional velocity optimization variables: observedSurfaceVelocityX, observedSurfaceVelocityY, observedSurfaceVelocityUncertainty, observedThicknessTendency, observedThicknessTendencyUncertainty, thicknessUncertainty')

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
if hasattr(fileout, 'history'):
   newhist = '\n'.join([thiscommand, getattr(fileout, 'history')])
else:
   newhist = thiscommand
setattr(fileout, 'history', newhist )

print("Completed creating land ice variables in new file. Now syncing to file.")
fileout.sync()

filein.close()
fileout.close()

print('\n** Successfully created {}.**'.format(options.fileoutName))
