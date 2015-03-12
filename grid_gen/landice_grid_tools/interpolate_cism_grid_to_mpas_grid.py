#!/usr/bin/env python
# Copy fields from a regular CISM grid into a pre-existing MPAS grid

import sys
import numpy as np
import netCDF4
from optparse import OptionParser
import math


print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n"
parser = OptionParser()
parser.description = "This script interpolates from a CISM grid to an MPAS grid using an ESMF weight interpolation file."
parser.add_option("-c", "--cism", dest="cismFile", help="CISM grid file to input.", default="cism.nc", metavar="FILENAME")
parser.add_option("-m", "--mpas", dest="mpasFile", help="MPAS grid file to output.", default="landice_grid.nc", metavar="FILENAME")
parser.add_option("-w", "--weight", dest="weightFile", help="ESMF weight file to input.  If not included, bilinear interpolation will be used", metavar="FILENAME")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

print "  CISM input file:  " + options.cismFile
print "  MPAS file to be modified:  " + options.mpasFile

if options.weightFile:
    interpType = 'weight'
    print "  Interpolation will be performed using ESMF weights from file:  " + options.weightFile
else:
    interpType = 'bilinear'
    print "  Bilinear interpolation will be performed."

print '' # make a space in stdout before further output


#----------------------------
# Define which time levels you want to use in the two files! (0-based indexing in python)
timelev = 0
timelevout = 0
#----------------------------

#----------------------------
# Map MPAS-CISM field names
fieldInfo = dict()
fieldInfo['thickness'] =     {'CISMname':'thk',  'scalefactor':1.0, 'CISMgrid':1, 'vertDim':False}
fieldInfo['bedTopography'] = {'CISMname':'topg', 'scalefactor':1.0, 'CISMgrid':1, 'vertDim':False}
fieldInfo['sfcMassBal'] =    {'CISMname':'acab', 'scalefactor':910.0/(3600.0*24.0*365.0), 'CISMgrid':1, 'vertDim':False}  # Assuming CISM density
fieldInfo['temperature'] =   {'CISMname':'tempstag', 'scalefactor':1.0, 'CISMgrid':1, 'vertDim':True}
fieldInfo['beta'] =          {'CISMname':'beta', 'scalefactor':1.0, 'CISMgrid':0, 'vertDim':False} # needs different mapping file...
#----------------------------

#----------------------------
#----------------------------
# Define needed functions
#----------------------------
#----------------------------

def ESMF_interp(sourceField):
    # Interpolates from the sourceField to the destinationField using ESMF weights
  try:
    # Initialize new field to 0 - required
    destinationField = np.zeros(xCell.shape)  # fields on cells only
    sourceFieldFlat = sourceField.flatten()  # Flatten source field
    for i in range(len(row)):
      destinationField[row[i]-1] = destinationField[row[i]-1] + S[i] * sourceFieldFlat[col[i]]
  except:
     'error in ESMF_interp'
  return destinationField

#----------------------------

def BilinearInterp(Value, CISMgridType):
    # Calculate bilinear interpolation of Value field from x, y to new ValueCell field (return value)  at xCell, yCell
    # This assumes that x, y, Value are regular CISM style grids and xCell, yCell, ValueCell are 1-D unstructured MPAS style grids
  try:
    if CISMgridType == 0:
        x = x0; y = y0
    elif CISMgridType == 1:
        x = x1; y = y1
    else:
        sys.exit('Error: unknown CISM grid type specified.')
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    ValueCell = np.zeros(xCell.shape)
    for i in range(len(xCell)):
       # Calculate the CISM grid cell indices (these are the lower index)
       xgrid = math.floor( (xCell[i]-x[0]) / dx )
       if xgrid >= len(x) - 1:
          xgrid = len(x) - 2
       elif xgrid < 0:
          xgrid = 0
       ygrid = math.floor( (yCell[i]-y[0]) / dy )
       if ygrid >= len(y) - 1:
          ygrid = len(y) - 2
       elif ygrid < 0:
          ygrid = 0
       #print xgrid, ygrid
       ValueCell[i] = Value[ygrid,xgrid] * (x[xgrid+1] - xCell[i]) * (y[ygrid+1] - yCell[i]) / (dx * dy) + \
                 Value[ygrid+1,xgrid] * (x[xgrid+1] - xCell[i]) * (yCell[i] - y[ygrid]) / (dx * dy) + \
                 Value[ygrid,xgrid+1] * (xCell[i] - x[xgrid]) * (y[ygrid+1] - yCell[i]) / (dx * dy) + \
                 Value[ygrid+1,xgrid+1] * (xCell[i] - x[xgrid]) * (yCell[i] - y[ygrid]) / (dx * dy) 
  except:
     'error in BilinearInterp'
  return ValueCell

#----------------------------

def interpolate_field(MPASfieldName):
    try:
        print '\n## %s ##'%MPASfieldName
        CISMfieldName = fieldInfo[MPASfieldName]['CISMname']
        CISMfield = CISMfile.variables[CISMfieldName][timelev,:,:]
        print '  CISM %s min/max:'%CISMfieldName, CISMfield.min(), CISMfield.max()
        # Call the appropriate routine for actually doing the interpolation
        if interpType == 'bilinear' or fieldInfo[MPASfieldName]['CISMgrid'] == 0:
            MPASfield = BilinearInterp(CISMfield, fieldInfo[MPASfieldName]['CISMgrid'])
        else:
            MPASfield = ESMF_interp(CISMfield)
        print '  interpolated MPAS %s min/max:'%MPASfieldName, MPASfield.min(), MPASfield.max()
#        # Don't let there be negative thickness
#        thickness[thickness<0.0] = 0.0
#        print 'cleaned MPAS thickness min/max:', thickness[:].min(), thickness[:].max()
        dims = MPASfile.variables[MPASfieldName].dimensions
        if 'Time' in dims:
            MPASfile.variables[MPASfieldName][timelevout,:] = MPASfield  # Time will always be leftmost index
        else:
            MPASfile.variables[MPASfieldName][:] = MPASfield
        if fieldInfo[MPASfieldName]['scalefactor'] != 1.0:
            MPASfield *= fieldInfo[MPASfieldName]['scalefactor']
            print '  scaled MPAS %s min/max:'%MPASfieldName, MPASfield.min(), MPASfield.max()

        del CISMfield, MPASfield
    except:
        print '  problem with %s field (e.g. not found in input file), skipping...'%CISMfieldName

#----------------------------
#----------------------------


#----------------------------
# Get weights from file
wfile = netCDF4.Dataset(options.weightFile, 'r')
S = wfile.variables['S'][:]
col = wfile.variables['col'][:]
row = wfile.variables['row'][:]
wfile.close()
#----------------------------



# Open the input file, get needed dimensions
try:
    CISMfile = netCDF4.Dataset(options.cismFile,'r')

    # Get the CISM dimensions if they exist
    try:
      level = len(CISMfile.dimensions['level'])
    except:
      print 'Input file is missing the dimension level.  Might not be a problem.'

    try:
      stagwbndlevel = len(CISMfile.dimensions['stagwbndlevel'])
    except:
      print 'Input file is missing the dimension stagwbndlevel.  Might not be a problem.'

    # Get CISM location variables if they exist
    try:
      x1 = CISMfile.variables['x1'][:]
      dx1 = x1[1] - x1[0]
      #print 'x1 min/max/dx:', x1.min(), x1.max(), dx1
      y1 = CISMfile.variables['y1'][:]
      dy1 = y1[1] - y1[0]
      #print 'y1 min/max/dx:', y1.min(), y1.max(), dy1

      ##x1 = x1 - (x1.max()-x1.min())/2.0  # This was for some shifted CISM grid but should not be used in general.
      ##y1 = y1 - (y1.max()-y1.min())/2.0
    except:
      print 'Input file is missing x1 and/or y1.  Might not be a problem.'
    
    try:
      x0 = CISMfile.variables['x0'][:]
      #print 'x0 min/max:', x0.min(), x0.max()
      y0 = CISMfile.variables['y0'][:]
      #print 'y0 min/max:', y0.min(), y0.max()

      ##x0 = x0 - (x0.max()-x0.min())/2.0
      ##y0 = y0 - (y0.max()-y0.min())/2.0

    except:
      print 'Input file is missing x0 and/or y0.  Might not be a problem.'

except:
    sys.exit('Error: The input file specified is either missing or lacking needed dimensions/variables.')



# Open the output file, get needed dimensions & variables
try:
    MPASfile = netCDF4.Dataset(options.mpasFile,'r+')
    try:
      nVertLevels = len(MPASfile.dimensions['nVertLevels'])
    except:
      print 'Output file is missing the dimension nVertLevels.  Might not be a problem.'

    try:
      # 1d vertical fields
      layerThicknessFractions = MPASfile.variables['layerThicknessFractions'][:]
    except:
      print 'Output file is missing the variable layerThicknessFractions.  Might not be a problem.'

    # '2d' spatial fields on cell centers
    xCell = MPASfile.variables['xCell'][:]
    #print 'xCell min/max:', xCell.min(), xCell.max()
    yCell = MPASfile.variables['yCell'][:]
    #print 'yCell min/max:', yCell.min(), yCell.max()

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
for field in fieldInfo:
    interpolate_field(field)


#try: 
#    inbeta = CISMfile.variables['beta']
#    beta = MPASfile.variables['betaTimeSeries']
#    print '\ninput beta min/max', inbeta[:].min(), inbeta[:].max()
#    #beta[timelevout,:] = BilinearInterp(x0, y0, inbeta[timelev,:,:], xCell, yCell)
#    beta[:,timelevout] = BilinearInterp(x0, y0, inbeta[timelev,:,:], xCell, yCell)
#    print 'interim beta min/max', beta[:].min(), beta[:].max()
#    # Make all beta be positive values
#    b = beta[:]
#    b[b<=0.0] = 1.0
#    beta[:] = b
#    print 'new beta min/max', beta[:].min(), beta[:].max()
#    del inbeta, beta, b
#except:
#    print '\nproblem with beta field (e.g. not found in input file), skipping...\n'


## If tempstag is present, then copy it directly over to the MPAS grid.  This currently requires that they have the same number of vertical levels.
#try: 
#    tempstag = CISMfile.variables['tempstag']
#    temperature = MPASfile.variables['temperature']
#    print '\ninput tempstag min/max', tempstag[:].min(), tempstag[:].max()
#    if stagwbndlevel == nVertLevels + 2:  # CISM includes the upper and lower boundaries as levels, MPAS does not.
#      #print range(1,level-1)
#      for i in range(1,stagwbndlevel-1):
#        print 'Copying level ', i+1, ' of ', stagwbndlevel , 'CISM staggered levels (ignoring CISM b.c. temp levels)'
#        temperature[timelevout,:,i-1] = BilinearInterp(x1, y1, tempstag[timelev,i,:,:], xCell, yCell)
#      print 'interim temperature min/max', temperature[:].min(), temperature[:].max()
#      # Don't let there be positive temperature
#      t = temperature[:]
#      t[t>0.0] = 0.0
#      temperature[:] = t
#      print 'new temperature min/max', temperature[:].min(), temperature[:].max()
#    else: 
#      raise Exception
#    del tempstag, temperature, t
#except:
#    print '\nproblem with tempstag field (e.g. not found in input file, differing number of layers), skipping...\n'


## If temp is present, copy over the level that is closest to the MPAS level.  Interpolation should be done, but starting with this simpler approach.  The CISM and MPAS grids need not have the same number of vertical levels.
#try: 
#    temp = CISMfile.variables['temp']
#    temperature = MPASfile.variables['temperature']
#    print '\ninput temp min/max', temp[:].min(), temp[:].max()
#    for i in range(nVertLevels):
#      print 'Copying level ', i+1, ' of ', nVertLevels
#      # grab the closest cism level (should interpolate)
#      if i==0:
#           icism=0
#      elif i==(nVertLevels-1):
#           icism=level-1
#      else:
#           icism=int(round( float(i+1)/float(nVertLevels) * float(level) ))-1
#      print 'Using CISM level ', icism+1
#      temperature[timelevout,:,i] = BilinearInterp(x1, y1, temp[timelev,icism,:,:], xCell, yCell)
#      print 'interim temperature min/max', temperature[:].min(), temperature[:].max()
#      # Don't let there be positive temperature
#      t = temperature[:]
#      t[t>0.0] = 0.0
#      temperature[:] = t
#      print 'new temperature min/max', temperature[:].min(), temperature[:].max()
#    else: 
#      raise Exception
#    del tempstag, temperature, t
#except:
#    print '\nproblem with tempstag field (e.g. not found in input file, differing number of layers), skipping...\n'



#try: 
#    balvel = CISMfile.variables['balvel']
#    observedSpeed = MPASfile.variables['observedSpeed']
#    print '\nbalvel min/max', balvel[:].min(), balvel[:].max()
#    observedSpeed[timelevout,:] = BilinearInterp(x0, y0, balvel[timelev,:,:], xCell, yCell)
#    print 'new observedSpeed min/max', observedSpeed[:].min(), observedSpeed[:].max()
#    del balvel, observedSpeed
#except:
#    print '\nproblem with balvel field (e.g. not found in input file), skipping...\n'
#

CISMfile.close()
MPASfile.close()

print '\nInterpolation completed.'
