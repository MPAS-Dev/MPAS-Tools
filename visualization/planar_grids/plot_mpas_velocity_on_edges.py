#!/usr/bin/env python
# Script capable of plotting most MPAS fields.  Invoke with 'python plot_mpas_field.py --help for details about how to use.
# Matt Hoffman, June 14, 2012

import sys, numpy, time
from optparse import OptionParser
import matplotlib.pyplot as plt
from netCDF4 import Dataset as NetCDFFile


print "** Gathering information.  (Invoke with 'python plot_mpas_field.py --help' for more details. All arguments are optional)"
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize; default: output.nc", metavar="FILE")
parser.add_option("-t", "--time", dest="time", help="time step to visualize (0 based); Give 'a' to animate all time levels or a single integer for a single plot of one time level; default: 0", metavar="TIME")
parser.add_option("-l", "--vertlevel", dest="vertlevel", help="vertical level to visualize (0 based); default: 0", metavar="LEVEL")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save animation as series of .png files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plot after animation is complete")
parser.add_option("-x", "--max", dest="maximum", help="maximum for color bar", metavar="MAX")
parser.add_option("-m", "--min", dest="minimum", help="minimum for color bar", metavar="MIN")
parser.add_option("-g", "--log", action="store_true", dest="log", help="include this flag to plot on log scale", metavar="LOG")
parser.add_option("-p", "--printval", action="store_true", dest="printval", help="include this flag to print numeric cell values (slow)", metavar="PRINTVAL")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

if not options.vertlevel:
	print "No vertical level provided. Using level 0 (surface), if level is needed."
        vert_level = 0
else:
        vert_level = int(options.vertlevel)
        print "Using vertical level " + options.vertlevel

f = NetCDFFile(options.filename,'r')

# Get grid stuff
try:
  xtime = f.variables['xtime'][:]  # Not needed unless trying to print actual time stamp
except:
  xtime = numpy.zeros((2,2))

if not options.time:
	print "No time provided. Using time 0, if time is needed."
        time_slices = 0
else:
        if options.time=='a':
           time_slices = numpy.arange( xtime.shape[0] )  # this is an easier way to get the time length than querying the dimension since the way to do that depends on which NetCDF module is being used...
           print "Animating all" + str(time_slices.size) + " time levels."
        else:
           time_slices = int(options.time)
           print "Using time level " +  options.time
try:
  xCell = f.variables['xCell']
  yCell = f.variables['yCell']
except:
  print 'xCell and/or yCell are missing or have problems.  Might not be a problem.'
try:
  xEdge = f.variables['xEdge']
  yEdge = f.variables['yEdge']
except:
  print 'xEdge and/or yEdge are missing or have problems.  Might not be a problem.'
try:
  angleEdge = f.variables['angleEdge']
except:
  print 'angleEdge is missing or has problems.  Might not be a problem.'
nCells = xCell.shape[0]  # This is an easier way to get nCells than querying the nCells dimension because the way to do that is module-specific


# This script needs certain things from the file.
varname = 'normalVelocity'
dynamicLimit = f.config_dynamic_thickness
thickness = f.variables['thickness']



# get the requested variable
if varname == 'uMag':
   uReconstructX = f.variables['uReconstructX']
   uReconstructY = f.variables['uReconstructY']
   var = numpy.zeros(uReconstructX.shape)
   var[:] = (uReconstructX[:] ** 2 + uReconstructY[:] **2)**0.5
   dims = uReconstructX.dimensions
else:
   var = f.variables[varname]
   dims = var.dimensions

# Determine what the appropriate X & Y values are:
if 'nCells' in dims:
   x = xCell[:]
   y = yCell[:]
elif 'nEdges' in dims:
   x = xEdge[:]
   y = yEdge[:]
elif 'nVertices' in dims:
   x = xVertex[:]
   y = yVertex[:]


# Get the needed slice and determine the plot title.  Make some assumptions about how dimensions are arranged:
plottitle = varname
if 'Time' in dims:
   time_length = var.shape[0]  # Assume time is the first dimension
   if 'nVertLevels' in dims:
      plottitle = plottitle + ', for layer ' + str(vert_level)
      var_slice = var[:,:,vert_level]
   else:
      var_slice = var[:,:]
else:
   print "Time is not a dimension of this variable.  Unable to animate it!"
   sys.exit()

if options.log:
        var_slice = numpy.log10(var_slice + 1.0e-5)
        plottitle = plottitle + ', log scale'

# Determine color axis max & min values to use.
if not options.maximum:
   maxval = var_slice[:].max()
else:
   maxval = float(options.maximum)

if not options.minimum:
   minval = var_slice[:].min()
else:
   minval = float(options.minimum)


# MAKE THE PLOT
print '** Beginning to create plot.'
plt.ion()
fig = plt.figure(1, facecolor='w', figsize=(7, 5), dpi=200)
ax = fig.add_subplot(111, aspect='equal')

# make an educated guess about how big the markers should be.
if nCells**0.5 < 60.0:
  markersize=  max( int(round(  3600.0/(nCells**0.5)  )), 1)
  markershape = 'h'  # use hexes if the points are big enough, otherwise just dots
else:
  markersize=  max( int(round(  1800.0/(nCells**0.5)  )), 1)
  markershape = '.'
print 'Using a markersize of ', markersize

if isinstance(time_slices, int) ==1:
   iterlist = [time_slices]
else:
   iterlist = time_slices

for t in iterlist:
    ax.cla()


    dynamicIndices = numpy.nonzero(thickness[:][t,:] >= dynamicLimit)[:]
    stagnantIndices = numpy.nonzero( numpy.logical_and( (thickness[:][t,:] < dynamicLimit), (thickness[:][t,:] > 0.0) ) )[:]

    plt.scatter(xCell[dynamicIndices], yCell[dynamicIndices], markersize, 'k', marker=markershape, edgecolors='none') # print dynamic ice locations with black hexagons
    plt.scatter(xCell[stagnantIndices], yCell[stagnantIndices], markersize, (0.5, 0.5, 0.5), marker=markershape, edgecolors='none') # print stagnant ice locations with gray hexagons

    nonzeroVelInd = numpy.nonzero(var_slice[t,:])
    plt.scatter(xEdge, yEdge, 3, 'k', marker='o', edgecolors='none') #print all vel's as block
    plt.scatter(xEdge[nonzeroVelInd], yEdge[nonzeroVelInd], 4, var_slice[t,nonzeroVelInd], marker='o', edgecolors='none', vmin=minval, vmax=maxval)  # print over the top nonzero vels with color

    if t==iterlist[0]:
       plt.colorbar()


    if options.printval:
       for iCell in range(x.shape[0]):
            plt.text(x[iCell], y[iCell], '{0:.1f}'.format(var_slice[t,iCell]), horizontalalignment='center', verticalalignment='center',  fontsize=3)
       dpi=300
    else:
       dpi=150

    plt.title( plottitle + ' at time ' + xtime[t,:].tostring().strip() ) #str(t) )
    plt.xlim( (x.min(), x.max() ) )
    plt.ylim( (y.min(), y.max() ) )
    plt.draw()
    time.sleep(0.05)
    if options.saveimages:
        plotname =  varname + '.' + '{0:04d}'.format(t) + '.'  + options.filename + '.png'
        #plt.ioff()
        plt.savefig(plotname,dpi=dpi)
        #plt.ion()
        print 'Saved plot as ' + plotname

plt.ioff()
print "Plotted " + varname + "."


if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     print 'Showing plot...  Close plot window to exit.'
     plt.show()

f.close()

