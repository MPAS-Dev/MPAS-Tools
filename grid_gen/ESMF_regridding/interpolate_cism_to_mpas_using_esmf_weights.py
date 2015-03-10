#!/usr/bin/env python
# Copy fields from a regular CISM grid into a pre-existing MPAS grid 

import sys, numpy
from netCDF4 import Dataset as NetCDFFile
import math
#import scipy.interpolate

#----------------------------
# Define which time levels you want to use in the two files! (0-based indexing in python)
timelev = 0
timelevout = 0
#----------------------------



# Check to see if the grid files were specified
if len(sys.argv) == 4:
  infilename = sys.argv[1]
  outfilename = sys.argv[2]
  weightfilename = sys.argv[3]
else:
  print '\nUsage:  python *.py [CISM_GRID.NC] [MPAS_GRID.NC] [weights.nc]'
  print len(sys.argv)
  sys.exit(0)


# get weights
wfile = NetCDFFile(weightfilename, 'r')
S = wfile.variables['S'][:]
col = wfile.variables['col'][:]
row = wfile.variables['row'][:]


# Open the input file, get needed dimensions
try:
    infile = NetCDFFile(infilename,'r')

    # Get the CISM dimensions if they exist
    try:
      level = len(infile.dimensions['level'])
    except:
      print 'Input file is missing the dimension level.  Might not be a problem.'

    try:
      stagwbndlevel = len(infile.dimensions['stagwbndlevel'])
    except:
      print 'Input file is missing the dimension stagwbndlevel.  Might not be a problem.'

#    try:
#      x1 = infile.dimensions['x1']
#      y1 = infile.dimensions['y1']
#    except:
#      print 'Input file is missing the dimensions x1, y1. '

#    try:
#      x1 = infile.dimensions['x1']
#      y1 = infile.dimensions['y1']
#    except:
#      print 'Input file is missing the dimensions x0, y0. '

    # Get CISM location variables if they exist
    try:
      x1 = infile.variables['x1'][:]
      dx1 = x1[1] - x1[0]
      print 'x1 min/max/dx:', x1.min(), x1.max(), dx1
      y1 = infile.variables['y1'][:]
      dy1 = y1[1] - y1[0]
      print 'y1 min/max/dx:', y1.min(), y1.max(), dy1

      ##x1 = x1 - (x1.max()-x1.min())/2.0  # This was for some shifted CISM grid but should not be used in general.
      ##y1 = y1 - (y1.max()-y1.min())/2.0
    except:
      print 'Input file is missing x1 and/or y1.  Might not be a problem.'
    
    try:
      x0 = infile.variables['x0'][:]
      print 'x0 min/max:', x0.min(), x0.max()
      y0 = infile.variables['y0'][:]
      print 'y0 min/max:', y0.min(), y0.max()

      ##x0 = x0 - (x0.max()-x0.min())/2.0
      ##y0 = y0 - (y0.max()-y0.min())/2.0

    except:
      print 'Input file is missing x0 and/or y0.  Might not be a problem.'

except:
    sys.exit('Error: The input file specified is either missing or lacking needed dimensions/variables.')



# Open the output file, get needed dimensions & variables
try:
    outfile = NetCDFFile(outfilename,'r+')
    try:
      nVertLevels = len(outfile.dimensions['nVertLevels'])
    except:
      print 'Output file is missing the dimension nVertLevels.  Might not be a problem.'

    try:
      # 1d vertical fields
      layerThicknessFractions = outfile.variables['layerThicknessFractions'][:]
    except:
      print 'Output file is missing the variable layerThicknessFractions.  Might not be a problem.'

    # '2d' spatial fields on cell centers
    xCell = outfile.variables['xCell'][:]
    print 'xCell min/max:', xCell.min(), xCell.max()
    yCell = outfile.variables['yCell'][:]
    print 'yCell min/max:', yCell.min(), yCell.max()

except:
    sys.exit('Error: The output grid file specified is either missing or lacking needed dimensions/variables.')


# Check the overlap of the grids
print '=================='
print 'CISM File extents:'
print '  x1 min, max:    ', x1.min(), x1.max()
print '  y1 min, max:    ', y1.min(), y1.max()
print 'MPAS File extents:'
print '  xCell min, max: ', xCell.min(), xCell.max()
print '  yCell min, max: ', yCell.min(), yCell.max()
print '=================='

#----------------------------
# try each field.  If it exists in the input file, it will be copied.  If not, it will be skipped.
# For now, include all variables defined for an MPAS land ice grid.  If more are added (e.g. SMB) they will need to be added below. 


thk = infile.variables['thk'][timelev,:,:]
print '\nthk min/max', thk[:].min(), thk[:].max()


thickness = numpy.zeros(xCell.shape)
thickness[:] = 0.0
for i in range(len(row)):
  thickness[row[i]-1] = thickness[row[i]-1] + S[i] * thk.flatten()[col[i]]
#    thickness[timelevout,:] = BilinearInterp(x1, y1, thk[timelev,:,:], xCell, yCell)
print 'interim thickness min/max', thickness[:].min(), thickness[:].max()

outfile.variables['thickness'][timelevout,:] = thickness

# Don't let there be negative thickness
#thickness[:] = thickness[:] * (thickness[:]>=0.0)
print 'new thickness min/max', thickness[:].min(), thickness[:].max()
del thk, thickness

outfile.close()

# do i=1, n_s
#   dst_field(row(i))=dst_field(row(i))+S(i)*src_field(col(i))
# enddo
