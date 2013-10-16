#!/usr/bin/python
# Script capable of plotting most MPAS fields on a planar mesh as a cross-section at a constant y location.  Invoke with 'python plot_mpas_field_xsect.py --help for details about how to use.
# Matt Hoffman, Jan. 2013

import numpy, time
from optparse import OptionParser
import matplotlib.pyplot as plt
from netCDF4 import Dataset as NetCDFFile


print "** Gathering information.  (Invoke with 'python plot_mpas_field_xsect.py --help' for more details. All arguments are optional)"
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize; default: output.nc", metavar="FILE")
parser.add_option("-v", "--var", dest="variable", help="variable to visualize; use 'uMag' for velocity magnitude; default: thickness", metavar="VAR")
parser.add_option("-y", "--y", dest="yval", help="y value to plot a cross-section on - script will find the closest valid value; use ncdump -v yCell to see your choices; leave empty for middle of domain", metavar="Y")
parser.add_option("-t", "--time", dest="time", help="time step to visualize (0 based); Give 'a' to animate all time levels or a single integer for a single plot of one time level; default: 0", metavar="TIME")
#parser.add_option("-l", "--level", dest="vertlevel", help="vertical level to visualize (0 based); default: 0", metavar="LEVEL")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plot as .png file")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")
parser.add_option("-x", "--max", dest="maximum", help="maximum for color bar", metavar="MAX")
parser.add_option("-m", "--min", dest="minimum", help="minimum for color bar", metavar="MIN")
parser.add_option("-g", "--log", action="store_true", dest="log", help="include this flag to plot on log scale", metavar="LOG")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

#if not options.vertlevel:
#	print "No vertical level provided. Using level 0 (surface), if level is needed."
#        vert_level = 0
#else:
#        vert_level = int(options.vertlevel)
#        print "Using vertical level " + options.vertlevel

if not options.variable:
        print "No variable provided. Using 'thickness'."
        varname = 'thickness'
else:
        varname = options.variable
        print "Using variable " + varname

f = NetCDFFile(options.filename,'r')

# Get grid stuff
xtime = f.variables['xtime'][:]  # Not needed unless trying to print actual time stamp

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
  print 'xCell and/or yCell is/are missing.  Might not be a problem.'
try:
  xEdge = f.variables['xEdge']
  yEdge = f.variables['yEdge']
except:
  print 'xEdge and/or yEdge is/are missing.  Might not be a problem.'
try:
  angleEdge = f.variables['angleEdge']
except:
  print 'angleEdge is missing.  Might not be a problem.'


# get the requested variable
if varname == 'uMag':
   uReconstructX = f.variables['uReconstructX']
   uReconstructY = f.variables['uReconstructY']
   var = numpy.zeros(uReconstructX.shape)
   var[:] = (uReconstructX[:] ** 2 + uReconstructY[:] **2)**0.5
   dims = uReconstructX.dimensions
elif varname == 'profile':
   var = f.variables['lowerSurface']
   dims = var.dimensions
   upperSurface = f.variables['upperSurface'][:]
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


if not options.yval:
        options.yval = (x.min() + x.max())/2.0
        print "No y level provided - using middle of domain: y=" + str(options.yval)
else:
        print "Using the closest y level to y=" + str(optiona.yval)


# Find the indices to the x-sect location
unique_ys=numpy.array(sorted(list(set(y))))
best_y=unique_ys[ numpy.absolute((unique_ys - float(options.yval))) == numpy.min(numpy.absolute(unique_ys - (float(options.yval)))) ]
print 'Found a best y value to use of:' + str(best_y)
theind=(y==best_y)

#def plotxsect():
#   """ this function actually plots the data"""


# Determine y-axis max & min values to use.
if not options.maximum:
   maxval = var[:].max()
   if varname == 'profile':
      maxval = max(maxval, upperSurface[:].max() )
else:
   maxval = float(options.maximum)

if not options.minimum:
   minval = var[:].min()
   if varname == 'profile':
      minval = min(minval, upperSurface[:].min() )
else:
   minval = float(options.minimum)

yaxisRng = maxval-minval


# MAKE THE PLOT
print '** Beginning to create plot.'
plt.ion()
fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111)

if isinstance(time_slices, int) ==1:
   iterlist = [time_slices]
else:
   iterlist = time_slices
#for t in numpy.arange(len(time_slices)):
for t in iterlist:
   ax.cla()

   # Get the needed slice and determine the plot title.  Make some assumptions about how dimensions are arranged:
   plottitle = varname + ' at y=' + str(best_y)
   if 'Time' in dims:
      plottitle = plottitle + ', at time ' + xtime[t,:].tostring().strip()
      var_slice = var[:][t,theind]
   else:
      var_slice = var[:][theind]
   if options.log:
        var_slice = numpy.log10(var_slice + 1.0e-5)
        plottitle = plottitle + ', log scale'

   # Make the plot
   plt.plot(xCell[:][theind], var_slice, '-o')
   if varname == 'profile':
     # also add the upper surface
     plt.plot(xCell[:][theind], upperSurface[t,theind], '-o')
   plt.title( plottitle )
   plt.ylim( ( minval - 0.1*yaxisRng, maxval + 0.1*yaxisRng) )
   plt.draw()
   print "Plotted " + varname + "."
   time.sleep(0.05)

   # save if needed.
   if options.saveimages:
        plotname =  varname + '.xsect.' + '{0:04d}'.format(t) + '.'  + options.filename + '.png' 

        plt.savefig(plotname)
        print 'Saved plot as ' + plotname

plt.ioff()

if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     print 'Showing plot...  Close plot window to exit.'
     plt.show()

if not isinstance(time_slices, int) and options.saveimages:
        print 'Try creating a movie with: ffmpeg -r 5 -i ' + varname + '.xsect.0%03d.' + options.filename + '.png ' + varname + '.xsect.' + options.filename + '.mp4'

f.close()

