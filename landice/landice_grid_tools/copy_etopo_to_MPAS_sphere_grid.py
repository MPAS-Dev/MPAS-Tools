#!/usr/bin/env python
# Copy fields from a regular CISM grid into a pre-existing MPAS grid 

import sys, numpy
from netCDF4 import Dataset
import math
from optparse import OptionParser

print "** Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.add_option("-e", "--etopo", dest="etopofile", help="etopo filename REQUIRED.  This currently has only been tested for the etopo5 dataset", metavar="FILENAME")
parser.add_option("-m", "--mpas", dest="mpasfile", help="mpas file to write to.  Defaults to 'landice_grid.nc'", metavar="FILENAME")
options, args = parser.parse_args()


if not options.etopofile:
    sys.exit('etopo filename is required.  Invoke with -h for more information.')
if not options.mpasfile:
    print "No MPAS filename specified, so using 'landice_grid.nc'."
    options.mpasfile = 'landice_grid.nc'


#----------------------------
# Define needed functions 
def BilinearInterp(x, y, Value, xCell, yCell):
    # Calculate bilinear interpolation of Value field from x, y to new ValueCell field (return value)  at xCell, yCell
    # This assumes that x, y, Value are regular CISM style grids and xCell, yCell, ValueCell are 1-D unstructured MPAS style grids
  try:  
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    ValueCell = numpy.zeros(xCell.shape)
    for i in range(len(xCell)):
       # Calculate the CISM grid cell indices (these are the lower index)
       xgrid = math.floor( (xCell[i]-x[0]) / dx )
       if xgrid >= len(x) - 1:
          xgrid = len(x) - 2
       ygrid = math.floor( (yCell[i]-y[0]) / dy )
       if ygrid >= len(y) - 1:
          ygrid = len(y) - 2
       #print xgrid, ygrid
       ValueCell[i] = Value[ygrid,xgrid] * (x[xgrid+1] - xCell[i]) * (y[ygrid+1] - yCell[i]) / (dx * dy) + \
                 Value[ygrid+1,xgrid] * (x[xgrid+1] - xCell[i]) * (yCell[i] - y[ygrid]) / (dx * dy) + \
                 Value[ygrid,xgrid+1] * (xCell[i] - x[xgrid]) * (y[ygrid+1] - yCell[i]) / (dx * dy) + \
                 Value[ygrid+1,xgrid+1] * (xCell[i] - x[xgrid]) * (yCell[i] - y[ygrid]) / (dx * dy) 
  except:
     'error in BilinearInterp'
  return  ValueCell
#----------------------------


# Open the input file, get variables
try:
    infile = Dataset(options.etopofile,'r')

    # Get the etopo dimensions if they exist
    elon = infile.variables['X']
    elat = infile.variables['Y']
    eZ = infile.variables['elev']
except:
    sys.exit('Error: The etopo input file specified is either missing or lacking needed dimensions/variables.')



# Open the output file, get needed dimensions & variables
try:
    outfile = Dataset(options.mpasfile, 'r+')
    # SHOULD CHECK TO MAKE SURE THIS IS ON A SPHERE!

    # '2d' spatial fields on cell centers
    latCell = outfile.variables['latCell']
    print 'latCell min/max:', latCell[:].min(), latCell[:].max()
    lonCell = outfile.variables['lonCell']
    print 'lonCell min/max:', lonCell[:].min(), lonCell[:].max()
    latCellD = latCell[:] * 180.0 / math.pi
    lonCellD = lonCell[:] * 180.0 / math.pi
except:
    sys.exit('Error: The output grid file specified is either missing or lacking needed dimensions/variables.')


#----------------------------
# try each field.  If it exists in the input file, it will be copied.  If not, it will be skipped.
# For now, include all variables defined for an MPAS land ice grid.  If more are added (e.g. SMB) they will need to be added below. 



if True:#try: 
    bedTopography = outfile.variables['bedTopography']
    print '\netopo elev min/max', eZ[:].min(), eZ[:].max()
    bedTopography[:] = BilinearInterp(elon[:], elat[:], eZ[:], lonCellD, latCellD)
    print 'new bedTopography min/max', bedTopography[:].min(), bedTopography[:].max()
#except:
#    print '\nproblem with elevation field (e.g. not found in input file), skipping...\n'
  

# Now smooth it!
coc = outfile.variables['cellsOnCell']
for iCell in range(len(bedTopography[:])):
   neighbors = coc[iCell,:]-1
   bedTopography[iCell] = bedTopography[:][neighbors].mean()



print '\nThis script is still experimental.  Make sure the in/out min/max values make sense before using the grid.'

infile.close()
outfile.close()


