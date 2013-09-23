#!/usr/bin/python
# Script to create a grid with land ice variables from an MPAS grid.
# I've only tested it with a periodic_hex grid, but it should work with any MPAS grid.
# Currently variable attributes are not copied (and periodic_hex does not assign any, so this is ok).  If variable attributes are added to periodic_hex, this script should be modified to copy them (looping over dir(var), skipping over variable function names "assignValue", "getValue", "typecode").

import sys, numpy
#try:
#  from Scientific.IO.NetCDF import NetCDFFile
#  netCDF_module = 'Scientific.IO.NetCDF'
#except ImportError:
#  try:
#    from netCDF4 import Dataset as NetCDFFile
#    netCDF_module = 'netCDF4'
#  except ImportError:
#      print 'Unable to import any of the following python modules:'
#      print '  Scientific.IO.NetCDF \n  netcdf4 '
#      print 'One of them must be installed.'
#      raise ImportError('No netCDF module found')

from netCDF4 import Dataset as NetCDFFile
netCDF_module = 'netCDF4'


# Check to see if a grid file was specified on the command line.
# If not, land_ice_grid.nc is used.
if len(sys.argv) > 1:
  if sys.argv[1][0] == '-': # The filename can't begin with a hyphen
    print '\nUsage:  python create_landice_grid_from_generic_MPAS_grid.py [GRID.NC]\nIf no filename is supplied, grid.nc will be used.'
    sys.exit(0)
  else:
    fileinName = sys.argv[1]
else:
  fileinName = 'grid.nc'



# Get some information about the input file
filein = NetCDFFile(fileinName,'r')
# vert_levs = filein.dimensions['nVertLevels']


# Define the new file to be output - this should perhaps be made a variant of the input file name
if netCDF_module == 'Scientific.IO.NetCDF':
    fileout = NetCDFFile("landice_grid.nc","w")
else:
    fileout = NetCDFFile("landice_grid.nc","w",format=filein.file_format)



# Copy over all the dimensions to the new file
# Note: looping over dimensions seems to result in them being written in seemingly random order.
#       I don't think this matters but it is not aesthetically pleasing.
#       It may be better to list them explicitly as I do for the grid variables, 
#       but this way ensures they all get included and is easier.
# Note: The UNLIMITED time dimension will return a dimension value of None with Scientific.IO.  This is what is supposed to happen.  See below for how to deal with assigning values to a variable with a unlimited dimension.  Special handling is needed with the netCDF module.
for dim in filein.dimensions.keys():
    if dim == 'nTracers': 
        pass  # Do nothing - we don't want this dimension 
    else:    # Copy over all other dimensions
      if netCDF_module == 'Scientific.IO.NetCDF':
        dimvalue = filein.dimensions[dim]
      else:
        if dim == 'Time':      
           dimvalue = None  # netCDF4 won't properly get this with the command below (you need to use the isunlimited method)
        else:
           dimvalue = len(filein.dimensions[dim])
      fileout.createDimension(dim, dimvalue)
# Create nVertLevelsPlus2 dimension
# fileout.createDimension('nVertLevelsPlus2', filein.dimensions['nVertLevels'] + 2)

# Create the dimensions needed for time-dependent forcings
# Note: These have been disabled in the fresh implementation of the landice core.  MH 9/19/13
#fileout.createDimension('nBetaTimeSlices', 1)
#fileout.createDimension('nSfcMassBalTimeSlices', 1)
#fileout.createDimension('nSfcAirTempTimeSlices', 1)
#fileout.createDimension('nBasalHeatFluxTimeSlices', 1)
#fileout.createDimension('nMarineBasalMassBalTimeSlices', 1)


# Copy over all of the required grid variables to the new file
vars2copy = ('latCell', 'lonCell', 'xCell', 'yCell', 'zCell', 'indexToCellID', 'latEdge', 'lonEdge', 'xEdge', 'yEdge', 'zEdge', 'indexToEdgeID', 'latVertex', 'lonVertex', 'xVertex', 'yVertex', 'zVertex', 'indexToVertexID', 'cellsOnEdge', 'nEdgesOnCell', 'nEdgesOnEdge', 'edgesOnCell', 'edgesOnEdge', 'weightsOnEdge', 'dvEdge', 'dcEdge', 'angleEdge', 'areaCell', 'areaTriangle', 'cellsOnCell', 'verticesOnCell', 'verticesOnEdge', 'edgesOnVertex', 'cellsOnVertex', 'kiteAreasOnVertex')
for varname in vars2copy:
   thevar = filein.variables[varname]
   if netCDF_module == 'Scientific.IO.NetCDF':
     datatype = thevar.typecode() 
   else:
     datatype = thevar.dtype
   newVar = fileout.createVariable(varname, datatype, thevar.dimensions)
   newVar[:] = thevar[:]
   # Create nVertLevelsPlus2 dimension - no longer used
       

# Create the land ice variables (all the shallow water vars in the input file can be ignored)
if netCDF_module == 'Scientific.IO.NetCDF':
    nVertLevels = fileout.dimensions['nVertLevels']
    datatype = 'd'
else:
    nVertLevels = len(filein.dimensions['nVertLevels'])
    datatype = filein.variables['xCell'].dtype  # Get the datatype for double precision float
#  Note: it may be necessary to make sure the Time dimension has size 1, rather than the 0 it defaults to.  For now, letting it be 0 which seems to be fine.
newvar = fileout.createVariable('layerThicknessFractions', datatype, ('nVertLevels', ))
newvar[:] = numpy.zeros(newvar.shape)
# Assign default values to layerThicknessFractions.  By default they will be uniform fractions.  Users can modify them in a subsequent step, but doing this here ensures the most likely values are already assigned. (Useful for e.g. setting up Greenland where the state variables are copied over but the grid variables are not modified.)
newvar[:] = 1.0 / nVertLevels


# With Scientific.IO.netCDF, entries are appended along the unlimited dimension one at a time by assigning to a slice.
# Therefore we need to assign to time level 0, and what we need to assign is a zeros array that is the shape of the new variable, exluding the time dimension!
newvar = fileout.createVariable('thickness', datatype, ('Time', 'nCells'))
newvar[0,:] = numpy.zeros( newvar.shape[1:] )   
newvar = fileout.createVariable('normalVelocity', datatype, ('Time', 'nEdges', 'nVertLevels'))
newvar[0,:,:] = numpy.zeros( newvar.shape[1:] )   
newvar = fileout.createVariable('temperature', datatype, ('Time', 'nCells', 'nVertLevels'))
newvar[0,:,:] = numpy.zeros( newvar.shape[1:] )   
# These landice variables are stored in the mesh currently, and therefore do not have a time dimension.
#    It may make sense to eventually move them to state.
newvar = fileout.createVariable('bedTopography', datatype, ('nCells',))
newvar[:] = numpy.zeros(newvar.shape)
newvar = fileout.createVariable('sfcMassBal', datatype, ('nCells',))
newvar[:] = numpy.zeros(newvar.shape)

# These boundary conditions are currently part of mesh, and are time independent.  If they change, make sure to adjust the dimensions here and in Registry.
# Note: These have been disabled in the fresh implementation of the landice core.  MH 9/19/13
#newvar = fileout.createVariable('betaTimeSeries', datatype, ( 'nCells', 'nBetaTimeSlices', ))
#newvar[:] = numpy.zeros(newvar.shape)
#newvar = fileout.createVariable('sfcMassBalTimeSeries', datatype, ( 'nCells', 'nSfcMassBalTimeSlices', ))
#newvar[:] = numpy.zeros(newvar.shape)
#newvar = fileout.createVariable('sfcAirTempTimeSeries', datatype, ( 'nCells', 'nSfcAirTempTimeSlices', ))
#newvar[:] = numpy.zeros(newvar.shape)
#newvar = fileout.createVariable('basalHeatFluxTimeSeries', datatype, ( 'nCells', 'nBasalHeatFluxTimeSlices',))
#newvar[:] = numpy.zeros(newvar.shape)
#newvar = fileout.createVariable('marineBasalMassBalTimeSeries', datatype, ( 'nCells', 'nMarineBasalMassBalTimeSlices',))
#newvar[:] = numpy.zeros(newvar.shape)

# Assign the global attributes
# Copy over the two attributes that are required by MPAS.  If any others exist in the input file, give a warning.
setattr(fileout, 'on_a_sphere', getattr(filein, 'on_a_sphere'))
setattr(fileout, 'sphere_radius', getattr(filein, 'sphere_radius'))
# If there are others that need to be copied, this script will need to be modified.  This Warning indicates that:
# Note: dir(file)  with Scientific.IO includes the following entries for functions in addition to the global attributes: 'close', 'createDimension', 'createVariable', 'flush', 'sync'.  NetCDF4 has a whole bunch more!
print "File had ", len(dir(filein)) - 5, "global attributes.  Copied on_a_sphere and sphere_radius."
print "Global attributes and functions: ", dir(filein)

filein.close()
fileout.close()

print 'Successfully created landice_grid.nc'
