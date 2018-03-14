#!/usr/bin/env python
'''
Script to plot common time-series from one or more landice globalStats files.
'''

import sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt



print "** Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser(description=__doc__)
parser.add_option("-1", dest="file1inName", help="input filename", default="globalStats.nc", metavar="FILENAME")
parser.add_option("-2", dest="file2inName", help="input filename", metavar="FILENAME")
parser.add_option("-3", dest="file3inName", help="input filename", metavar="FILENAME")
parser.add_option("-4", dest="file4inName", help="input filename", metavar="FILENAME")
options, args = parser.parse_args()


# create axes to plot into
fig = plt.figure(1, facecolor='w')

nrow=4
ncol=2

#xtickSpacing = 20.0

axVol = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
plt.ylabel('volume (m$^3$)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()
axX = axVol

axVAF = fig.add_subplot(nrow, ncol, 2, sharex=axX)
plt.xlabel('Year')
plt.ylabel('VAF (m^3)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axVolGround = fig.add_subplot(nrow, ncol, 3, sharex=axX)
plt.xlabel('Year')
plt.ylabel('grounded volume (m$^3$)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axVolFloat = fig.add_subplot(nrow, ncol, 4, sharex=axX)
plt.xlabel('Year')
plt.ylabel('floating volume (m$^3$)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axGrdArea = fig.add_subplot(nrow, ncol, 5, sharex=axX)
plt.xlabel('Year')
plt.ylabel('grounded area (m$^2$)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axFltArea = fig.add_subplot(nrow, ncol, 6, sharex=axX)
plt.xlabel('Year')
plt.ylabel('floating area (m$^2$)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axGLflux = fig.add_subplot(nrow, ncol, 7, sharex=axX)
plt.xlabel('Year')
plt.ylabel('GL flux (kg/yr)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axCalvFlux = fig.add_subplot(nrow, ncol, 8, sharex=axX)
plt.xlabel('Year')
plt.ylabel('calving flux (kg/yr)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()



def plotStat(fname):
    print "Reading and plotting file: " + fname

    name = fname

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0

    vol = f.variables['totalIceVolume'][:]
    axVol.plot(yr, vol, label=name)

    VAF = f.variables['volumeAboveFloatation'][:]
    axVAF.plot(yr, VAF, label=name)

    volGround = f.variables['groundedIceVolume'][:]
    axVolGround.plot(yr, volGround, label=name)

    volFloat = f.variables['floatingIceVolume'][:]
    axVolFloat.plot(yr, volFloat, label=name)

    areaGrd = f.variables['groundedIceArea'][:]
    axGrdArea.plot(yr, areaGrd, label=name)

    areaFlt = f.variables['floatingIceArea'][:]
    axFltArea.plot(yr, areaFlt, label=name)

    GLflux = f.variables['groundingLineFlux'][:]
    axGLflux.plot(yr, GLflux, label=name)

    calvFlux = f.variables['totalCalvingFlux'][:]
    axCalvFlux.plot(yr, calvFlux, label=name)


    f.close()


plotStat(options.file1inName)


if(options.file2inName):
   plotStat(options.file2inName)

if(options.file3inName):
   plotStat(options.file3inName)

if(options.file4inName):
   plotStat(options.file4inName)

axCalvFlux.legend(loc='best', prop={'size': 6})

print "Generating plot."
plt.show()


