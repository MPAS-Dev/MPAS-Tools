#!/usr/bin/env python
'''
Script to plot common time-series from one or more landice globalStats files.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt

rhoi = 910.0


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-1", dest="file1inName", help="input filename", default="globalStats.nc", metavar="FILENAME")
parser.add_option("-2", dest="file2inName", help="input filename", metavar="FILENAME")
parser.add_option("-3", dest="file3inName", help="input filename", metavar="FILENAME")
parser.add_option("-4", dest="file4inName", help="input filename", metavar="FILENAME")
parser.add_option("-5", dest="file5inName", help="input filename", metavar="FILENAME")
parser.add_option("-6", dest="file6inName", help="input filename", metavar="FILENAME")
parser.add_option("-7", dest="file7inName", help="input filename", metavar="FILENAME")
parser.add_option("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="Gt", metavar="FILENAME")
parser.add_option("-c", dest="plotChange", help="plot time series as change from initial.  (not applied to GL flux or calving flux)  Without this option, the full magnitude of time series is used", action='store_true')
options, args = parser.parse_args()

print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))

# create axes to plot into
fig = plt.figure(1, figsize=(9, 11), facecolor='w')

nrow=4
ncol=2

#xtickSpacing = 20.0

if options.units == "m3":
   massUnit = "m$^3$"
elif options.units == "kg":
   massUnit = "kg"
elif options.units == "Gt":
   massUnit = "Gt"
else:
   sys.exit("Unknown mass/volume units")
print("Using volume/mass units of: ", massUnit)

axVol = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
if options.plotChange:
   plt.ylabel('volume change ({})'.format(massUnit))
else:
   plt.ylabel('volume ({})'.format(massUnit))
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()
axX = axVol

axVAF = fig.add_subplot(nrow, ncol, 2, sharex=axX)
plt.xlabel('Year')
if options.plotChange:
   plt.ylabel('VAF change ({})'.format(massUnit))
else:
   plt.ylabel('VAF ({})'.format(massUnit))
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axVolGround = fig.add_subplot(nrow, ncol, 3, sharex=axX)
plt.xlabel('Year')
if options.plotChange:
   plt.ylabel('grounded volume change ({})'.format(massUnit))
else:
   plt.ylabel('grounded volume ({})'.format(massUnit))
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axVolFloat = fig.add_subplot(nrow, ncol, 4, sharex=axX)
plt.xlabel('Year')
if options.plotChange:
   plt.ylabel('floating volume change ({})'.format(massUnit))
else:
   plt.ylabel('floating volume ({})'.format(massUnit))
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axGrdArea = fig.add_subplot(nrow, ncol, 5, sharex=axX)
plt.xlabel('Year')
if options.plotChange:
   plt.ylabel('grounded area change (m$^2$)')
else:
   plt.ylabel('grounded area (m$^2$)')
#plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axFltArea = fig.add_subplot(nrow, ncol, 6, sharex=axX)
plt.xlabel('Year')
if options.plotChange:
   plt.ylabel('floating area change (m$^2$)')
else:
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
    print("Reading and plotting file: {}".format(fname))

    name = fname

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    yr = yr-yr[0]
    dt = f.variables['deltat'][:]/3.15e7
    print(yr.max())

    vol = f.variables['totalIceVolume'][:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       vol = vol * rhoi
    elif options.units == "Gt":
       vol = vol * rhoi / 1.0e12
    if options.plotChange:
        vol = vol - vol[0]
    axVol.plot(yr, vol, label=name)

    VAF = f.variables['volumeAboveFloatation'][:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       VAF = VAF * rhoi
    elif options.units == "Gt":
       VAF = VAF * rhoi / 1.0e12
    if options.plotChange:
        VAF = VAF - VAF[0]
    axVAF.plot(yr, VAF, label=name)

    volGround = f.variables['groundedIceVolume'][:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       volGround = volGround * rhoi
    elif options.units == "Gt":
       volGround = volGround * rhoi / 1.0e12
    if options.plotChange:
        volGround = volGround - volGround[0]
    axVolGround.plot(yr, volGround, label=name)

    volFloat = f.variables['floatingIceVolume'][:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       volFloat = volFloat * rhoi
    elif options.units == "Gt":
       volFloat = volFloat * rhoi / 1.0e12
    if options.plotChange:
        volFloat = volFloat - volFloat[0]
    axVolFloat.plot(yr, volFloat, label=name)

    areaGrd = f.variables['groundedIceArea'][:]
    if options.plotChange:
        areaGrd = areaGrd - areaGrd[0]
    axGrdArea.plot(yr, areaGrd, label=name)

    areaFlt = f.variables['floatingIceArea'][:]
    if options.plotChange:
        areaFlt = areaFlt - areaFlt[0]
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

if(options.file5inName):
   plotStat(options.file5inName)

if(options.file6inName):
   plotStat(options.file6inName)

if(options.file7inName):
   plotStat(options.file7inName)

axCalvFlux.legend(loc='best', prop={'size': 6})

print("Generating plot.")
fig.tight_layout()
plt.show()

