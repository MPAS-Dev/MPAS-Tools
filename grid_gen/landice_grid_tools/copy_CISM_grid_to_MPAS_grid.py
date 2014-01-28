#!/usr/bin/python
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
if len(sys.argv) == 3:
  infilename = sys.argv[1]
  outfilename = sys.argv[2]
else:
  print '\nUsage:  python *.py [CISM_GRID.NC] [MPAS_GRID.NC]'
  print len(sys.argv)
  sys.exit(0)

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
    except:
      print 'Input file is missing x1 and/or y1.  Might not be a problem.'
    
    try:
      x0 = infile.variables['x0'][:]
      print 'x0 min/max:', x0.min(), x0.max()
      y0 = infile.variables['y0'][:]
      print 'y0 min/max:', y0.min(), y0.max()
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

try: 
    thk = infile.variables['thk']
    thickness = outfile.variables['thickness']
    print '\nthk min/max', thk[:].min(), thk[:].max()
    thickness[timelevout,:] = BilinearInterp(x1, y1, thk[timelev,:,:], xCell, yCell)
    print 'interim thickness min/max', thickness[:].min(), thickness[:].max()
    # Don't let there be negative thickness
    thickness[:] = thickness[:] * (thickness[:]>=0.0)
    print 'new thickness min/max', thickness[:].min(), thickness[:].max()
    del thk, thickness
except:
    print '\nproblem with thk field (e.g. not found in input file), skipping...\n'


try: 
    topg = infile.variables['topg']
    bedTopography = outfile.variables['bedTopography']
    print '\ntopg min/max', topg[:].min(), topg[:].max()
    bedTopography[:] = BilinearInterp(x1, y1, topg[timelev,:,:], xCell, yCell)
    print 'new bedTopography min/max', bedTopography[:].min(), bedTopography[:].max()
    del topg, bedTopography
except:
    print '\nproblem with topg field (e.g. not found in input file), skipping...\n'
  

try: 
    inbeta = infile.variables['beta']
    beta = outfile.variables['betaTimeSeries']
    print '\ninput beta min/max', inbeta[:].min(), inbeta[:].max()
    #beta[timelevout,:] = BilinearInterp(x0, y0, inbeta[timelev,:,:], xCell, yCell)
    beta[:,timelevout] = BilinearInterp(x0, y0, inbeta[timelev,:,:], xCell, yCell)
    print 'interim beta min/max', beta[:].min(), beta[:].max()
    # Make all beta be positive values
    b = beta[:]
    b[b<=0.0] = 1.0
    beta[:] = b
    print 'new beta min/max', beta[:].min(), beta[:].max()
    del inbeta, beta, b
except:
    print '\nproblem with beta field (e.g. not found in input file), skipping...\n'


# If tempstag is present, then copy it directly over to the MPAS grid.  This currently requires that they have the same number of vertical levels.
try: 
    tempstag = infile.variables['tempstag']
    temperature = outfile.variables['temperature']
    print '\ninput tempstag min/max', tempstag[:].min(), tempstag[:].max()
    if stagwbndlevel == nVertLevels + 2:  # CISM includes the upper and lower boundaries as levels, MPAS does not.
      #print range(1,level-1)
      for i in range(1,stagwbndlevel-1):
        print 'Copying level ', i+1, ' of ', stagwbndlevel , 'CISM staggered levels (ignoring CISM b.c. temp levels)'
        temperature[timelevout,:,i-1] = BilinearInterp(x1, y1, tempstag[timelev,i,:,:], xCell, yCell)
      print 'interim temperature min/max', temperature[:].min(), temperature[:].max()
      # Don't let there be positive temperature
      t = temperature[:]
      t[t>0.0] = 0.0
      temperature[:] = t
      print 'new temperature min/max', temperature[:].min(), temperature[:].max()
    else: 
      raise Exception
    del tempstag, temperature, t
except:
    print '\nproblem with tempstag field (e.g. not found in input file, differing number of layers), skipping...\n'


# If temp is present, copy over the level that is closest to the MPAS level.  Interpolation should be done, but starting with this simpler approach.  The CISM and MPAS grids need not have the same number of vertical levels.
try: 
    temp = infile.variables['temp']
    temperature = outfile.variables['temperature']
    print '\ninput temp min/max', temp[:].min(), temp[:].max()
    for i in range(nVertLevels):
      print 'Copying level ', i+1, ' of ', nVertLevels
      # grab the closest cism level (should interpolate)
      if i==0:
           icism=0
      elif i==(nVertLevels-1):
           icism=level-1
      else:
           icism=int(round( float(i+1)/float(nVertLevels) * float(level) ))-1
      print 'Using CISM level ', icism+1
      temperature[timelevout,:,i] = BilinearInterp(x1, y1, temp[timelev,icism,:,:], xCell, yCell)
      print 'interim temperature min/max', temperature[:].min(), temperature[:].max()
      # Don't let there be positive temperature
      t = temperature[:]
      t[t>0.0] = 0.0
      temperature[:] = t
      print 'new temperature min/max', temperature[:].min(), temperature[:].max()
    else: 
      raise Exception
    del tempstag, temperature, t
except:
    print '\nproblem with tempstag field (e.g. not found in input file, differing number of layers), skipping...\n'



try: 
    balvel = infile.variables['balvel']
    observedSpeed = outfile.variables['observedSpeed']
    print '\nbalvel min/max', balvel[:].min(), balvel[:].max()
    observedSpeed[timelevout,:] = BilinearInterp(x0, y0, balvel[timelev,:,:], xCell, yCell)
    print 'new observedSpeed min/max', observedSpeed[:].min(), observedSpeed[:].max()
    del balvel, observedSpeed
except:
    print '\nproblem with balvel field (e.g. not found in input file), skipping...\n'
  

try: 
    acab = infile.variables['acab']
    smb = outfile.variables['sfcMassBalTimeSeries']
    print '\nacab min/max', acab[:].min(), acab[:].max()
    smbmpas = BilinearInterp(x1, y1, acab[timelev,:,:], xCell, yCell)
    for t in range(smb.shape[1]):
            print 'Inserting CISM acab from time ' + str(timelev) + ' into MPAS sfcMassBalTimeSeries time index ' + str(t)
            smb[:,t] = smbmpas[:]
    print 'new smb min/max', smb[:].min(), smb[:].max()
    del acab, smb
except:
    print '\nproblem with acab field (e.g. not found in input file), skipping...\n'



# Variables on edges...?
#normalVelocity = outfile.variables['normalVelocity'][:]


# Old implementation, including the super-slow scipy.interpolate object
#try:
#    thk = infile.variables['thk'][:]
#    print 'thk min/max', thk.min(), thk.max()
#    #f = scipy.interpolate.RectBivariateSpline(y1, x1, thk[timelev, :, :], kx=1, ky=1)
#    #print 'got f'
#    #znew = f(yCell, xCell)
#    for i in range(len(xCell)):
#       # Calculate the CISM grid cell indices (these are the lower index)
#       xgrid = math.floor( xCell[i] / dx1 )
#       if xgrid >= len(x1)-1:
#          xgrid = len(x1) - 2
#       ygrid = math.floor( yCell[i] / dy1 )
#       if ygrid >= len(y1)-1:
#          ygrid = len(y1) - 2
#       #print xgrid, ygrid
#       thickness[0,i] = thk[timelev,ygrid,xgrid] * (x1[xgrid+1] - xCell[i]) * (y1[ygrid+1] - yCell[i]) / (dx1 * dy1) + \
#                 thk[timelev,ygrid+1,xgrid] * (x1[xgrid+1] - xCell[i]) * (yCell[i] - y1[ygrid]) / (dx1 * dy1) + \
#                 thk[timelev,ygrid,xgrid+1] * (xCell[i] - x1[xgrid]) * (y1[ygrid+1] - yCell[i]) / (dx1 * dy1) + \
#                 thk[timelev,ygrid+1,xgrid+1] * (xCell[i] - x1[xgrid]) * (yCell[i] - y1[ygrid]) / (dx1 * dy1) 
#    print 'new thickness min/max', thickness[:].min(), thickness[:].max()
#except:
#    print 'problem with thk field (e.g. not found in input file), skipping...'

print '\nThis script is still experimental.  Make sure the in/out min/max values make sense before using the grid.'

infile.close()
outfile.close()


