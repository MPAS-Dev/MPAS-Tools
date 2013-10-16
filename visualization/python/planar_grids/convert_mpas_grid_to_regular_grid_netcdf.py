#!/usr/bin/python
# Script to re-grid variables in a planar MPAS output file to 
# a regular rectangular grid in a new file.  
# Currently, the script is setup to regrid any 2d and 3d variables that have dimensions of both Time
# and nCells.  The output filename has the '.regulargrid.nc' appended to it.  
# It uses a Natural Neighbors algorithm.
# All global and variable attributes are copied over.  Supporing Scientific.IO.NetCDF requires doing this in a clumsy way.  The netCDF4 module allows the use of the .ncattrs() method to get just the NetCDF attributes and none of the python attributes for those objects, but this is not implemented for cross-compatibility.
# Matt Hoffman, LANL, June 28, 2012

try:
  from Scientific.IO.NetCDF import NetCDFFile
  netCDF_module = 'Scientific.IO.NetCDF'
except ImportError:
  try:
    from netCDF4 import Dataset as NetCDFFile
    netCDF_module = 'netCDF4'
  except ImportError:
      print 'Unable to import any of the following python modules:'
      print '  Scientific.IO.NetCDF \n  netcdf4 '
      print 'One of them must be installed.'
      raise ImportError('No netCDF module found')
from optparse import OptionParser
import numpy, matplotlib.delaunay

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to re-grid", metavar="FILE")
parser.add_option("-x", "--nx", dest="nx", help="number of output cells to use in x direction (optional)", metavar="NX")
parser.add_option("-y", "--ny", dest="ny", help="number of output cells to use in x direction (optional)", metavar="NY")
parser.add_option("-r", "--range", dest="range", help="comma delimited list of min,max x,y values to use for the output grid (optional).  Format: '-r xmin,xmax,ymin,ymax'", metavar="RANGE")

options, args = parser.parse_args()
if not options.filename:
	parser.error("Filename is a required input.")

print "This script will attemp to re-grid all variables in the input file that contain both dimensions nCells and Time to a regular grid that can be viewed in, e.g, ncview.  Invoke with 'python convert_mpas_grid_to_regular_grid_netcdf.py --help' for more information. \n"

# Get some information about the input file
filein = NetCDFFile(options.filename,'r')
xIn = filein.variables['xCell'][:]
yIn = filein.variables['yCell'][:]
times = filein.variables['xtime']
time_length = times.shape[0]
if netCDF_module == 'Scientific.IO.NetCDF':
  vert_levs = filein.dimensions['nVertLevels']
else:
  vert_levs = len(filein.dimensions['nVertLevels'])

# Determine what to use for the grid min/max values
#if options.range:
try:
  xmin, xmax, ymin, ymax = map(float, options.range.split(','))
#else:
except:
  xmin = xIn.min()
  xmax = xIn.max()
  ymin = yIn.min()
  ymax = yIn.max()

# Make a guess at the proper dimensions of the output file 
dcedge = filein.variables['dcEdge']
resolution = dcedge[:].max()
# Compute nx and ny from the info in the input file, or use values supplied on command line
if options.nx:
  nx = int(options.nx)
else:
  nx = int((xmax - xmin)/resolution - 1)
if options.ny:
  ny = int(options.ny)
else:
  ny = int((ymax - ymin)/resolution)


# ==== Define the new file to be output - this should be made a variant of the input file name
fileoutname = options.filename + ".regulargrid.nc"
fileout = NetCDFFile(fileoutname,"w")

# ====Create dimensions in output file
# Create the new x,y dimensions for the new file
fileout.createDimension('x',nx)
fileout.createDimension('y',ny)
# Copy over all the dimensions to the new file
for dim in filein.dimensions.keys():
    # print 'DIMENSION: ', dim
    # print 'HAS VALUE: ', filein.dimensions[dim]
    if netCDF_module == 'Scientific.IO.NetCDF':
      fileout.createDimension(dim, filein.dimensions[dim])
    else:
      fileout.createDimension(dim, len(filein.dimensions[dim]))


# ====Copy over global attributes
for a in dir(filein):
  if not ( any(x in a for x in ('close', 'createDimension', 'createVariable', 'flush', 'sync', 'groups', 'dimensions', 'variables', 'dtype', 'file_format', '_nunlimdim', 'path', 'parent', 'ndim', 'maskandscale', 'cmptypes', 'vltypes', 'createCompoundType', 'createGroup', 'createVLType', 'delncattr', 'getncattr', 'maskanscale', 'ncattrs', 'renameDimension', 'renameVariable', 'setncattr') ) or ('_' in a) ):  # don't copy these
     setattr(fileout, a, getattr(filein, a) )

# ====Create x,y,time coordinate variables
try: 
  print 'Attempting to copy time coordinate variable'
  newTime = fileout.createVariable('xtime', 'c', ('Time', 'StrLen') )
  # Copy over attributes
  for a in dir(filein.variables['xtime']):
     if not (any(x in a for x in ('typecode', 'assignValue', 'chunking', 'delncattr', 'dimensions', 'dtype', 'endian', 'filters', 'getValue', 'get_var_chunk_cache', 'getncattr', 'group', 'maskandscale', 'ncattrs', 'ndim', 'set_auto_maskandscale', 'set_var_chunk_cache', 'setncattr', 'shape', 'size') ) or ('_' in a) ):  # don't copy these
         setattr(newTime, a, getattr(filein.variables['xtime'], a) )
except:
  print 'Error in copying time field.  Skipping it...'

print 'Attempting to create x, y coordinate variables'
newX = fileout.createVariable('x', 'd', ('x',) )
newX[:] = numpy.linspace(xmin, xmax, nx)
newY = fileout.createVariable('y', 'd', ('y',) )
newY[:] = numpy.linspace(ymin, ymax, ny)

# ====Loop over all variables in file and find the ones we want to convert
# Create a triang object that can be used for all variables (to speed the conversion process)
triang = matplotlib.delaunay.Triangulation(xIn, yIn)
for var in filein.variables.keys():
   # print var
   thevar = filein.variables[var]
   dimflag = 0
   # Look for variables that have these two dimensions
   if netCDF_module == 'Scientific.IO.NetCDF':
     typecode = thevar.typecode()
   else:
     typecode = thevar.dtype
   if (  all(d in thevar.dimensions for d in ('nCells', 'Time') )  and  any(tc in typecode for tc in ('d', 'f') )  ):
       print 'Attempting to re-grid variable: ', var
       # Determine proper dimensions for the new variable
       dims2use = ()
       lastdim = 0
       for dim in thevar.dimensions:
           if dim == 'nCells':
              dims2use = dims2use + ('y', 'x')
           else:
              dims2use = dims2use + (dim,)
           if dim == 'nVertLevels':
              lastdim = vert_levs
           elif dim == 'nVertLevelsPlus2':
              lastdim = vert_levs + 2
       # Create variable
       newVar = fileout.createVariable(var,typecode,dims2use)
       # Copy over attributes
       for a in dir(thevar):
          if not (any(x in a for x in ('typecode', 'assignValue', 'chunking', 'delncattr', 'dimensions', 'dtype', 'endian', 'filters', 'getValue', 'get_var_chunk_cache', 'getncattr', 'group', 'maskandscale', 'ncattrs', 'ndim', 'set_auto_maskandscale', 'set_var_chunk_cache', 'setncattr', 'shape', 'size') ) or ('_' in a) ):  # don't copy these
             setattr(newVar, a, getattr(thevar, a) )
             
       # Loop over time, and vert levels if needed, interpolating each 2d field for this variable
       for t in range(0,time_length):      
          if 'nVertLevels' in thevar.dimensions:
             for v in range(0,vert_levs):
                originalvalue = thevar[t,:,v]
                # Create an interpolator object for this variable
                interpolator = matplotlib.delaunay.NNInterpolator(triang, originalvalue, default_value=0.0)
                # use the object to regrid
                varGridded = interpolator[ymin:ymax:ny*1j, xmin:xmax:nx*1j]
                newVar[t,:,:,v] = varGridded[:,:]
          else:
             originalvalue = thevar[t,:]
             # Create an interpolator object for this variable
             interpolator = matplotlib.delaunay.NNInterpolator(triang, originalvalue, default_value=0.0)
             # use the object to regrid
             varGridded = interpolator[ymin:ymax:ny*1j, xmin:xmax:nx*1j]            
             newVar[t,:,:] = varGridded[:,:]

print 'Saved re-gridded fields to: ', fileoutname
 

fileout.sync()
filein.close()
fileout.close()
