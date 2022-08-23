#!/usr/bin/env python
'''
Script to plot common time-series from one or more landice regionalStats files.
Matt Hoffman, 8/23/2022
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
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
parser.add_option("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="Gt", metavar="FILENAME")
parser.add_option("-n", dest="fileRegionNames", help="region name filename", metavar="FILENAME")
options, args = parser.parse_args()

print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))


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

if options.fileRegionNames:
   fn = Dataset(options.fileRegionNames, 'r')
   rNames = fn.variables['regionNames'][:]

# Get nRegions from first file
f = Dataset(options.file1inName, 'r')
nRegions = len(f.dimensions['nRegions'])

nrow=4
ncol=4
if nRegions > nrow*ncol:
    sys.exit("ERROR: Number of regions exceeds number of plots.  Please adjust nrow, ncol as needed.")

# Set up Figure 1: volume stats overview
fig1, axs1 = plt.subplots(nrow, ncol, figsize=(13, 11), num=1)
fig1.suptitle('Mass change summary', fontsize=14)
for reg in range(nRegions):
   plt.sca(axs1.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('volume change ({})'.format(massUnit))
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs1.flatten()[reg].set_title(rNames[reg,:].tobytes().decode('utf-8').strip())
   if reg == 0:
      axX = axs1.flatten()[reg]
   else:
      axs1.flatten()[reg].sharex(axX)

# Set up Figure 2: grd MB
fig2, axs2 = plt.subplots(nrow, ncol, figsize=(13, 11), num=2)
fig2.suptitle('Grounded mass change', fontsize=14)
for reg in range(nRegions):
   plt.sca(axs2.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('volume change ({})'.format(massUnit))
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs2.flatten()[reg].set_title(rNames[reg,:].tobytes().decode('utf-8').strip())
   if reg == 0:
      axX = axs2.flatten()[reg]
   else:
      axs2.flatten()[reg].sharex(axX)

# Set up Figure 3: flt MB
fig3, axs3 = plt.subplots(nrow, ncol, figsize=(13, 11), num=3)
fig3.suptitle('Floating mass change', fontsize=14)
for reg in range(nRegions):
   plt.sca(axs3.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('volume change ({})'.format(massUnit))
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs3.flatten()[reg].set_title(rNames[reg,:].tobytes().decode('utf-8').strip())
   if reg == 0:
      axX = axs3.flatten()[reg]
   else:
      axs3.flatten()[reg].sharex(axX)




def plotStat(fname, sty):
    print("Reading and plotting file: {}".format(fname))

    name = fname

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    nRegionsLocal = len(f.dimensions['nRegions'])
    if nRegionsLocal != nRegions:
        sys.exit(f"ERROR: Number of regions in file {fname} does not match number of regions in first input file!")

    # Fig 1: summary plot
    vol = f.variables['regionalIceVolume'][:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       vol = vol * rhoi
    elif options.units == "Gt":
       vol = vol * rhoi / 1.0e12
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, vol[:,r] - vol[0,r], label="total", linestyle=sty, color='k')

    VAF = f.variables['regionalVolumeAboveFloatation'][:]
    VAF = VAF[:,:] - VAF[0,:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       VAF = VAF * rhoi
    elif options.units == "Gt":
       VAF = VAF * rhoi / 1.0e12
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, VAF[:,r] - VAF[0,r], label="VAF", linestyle=sty, color='m')

    volGround = f.variables['regionalGroundedIceVolume'][:]
    volGround = volGround[:,:] - volGround[0,:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       volGround = volGround * rhoi
    elif options.units == "Gt":
       volGround = volGround * rhoi / 1.0e12
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, volGround[:,r] - volGround[0,r], label="grd", linestyle=sty, color='b')

    volFloat = f.variables['regionalFloatingIceVolume'][:]
    volFloat = volFloat[:,:] - volFloat[0,:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       volFloat = volFloat * rhoi
    elif options.units == "Gt":
       volFloat = volFloat * rhoi / 1.0e12
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, volFloat[:,r] - volFloat[0,r], label="flt", linestyle=sty, color='g')

#    areaGrd = f.variables['groundedIceArea'][:]
#    axGrdArea.plot(yr, areaGrd, label=name)
#
#    areaFlt = f.variables['floatingIceArea'][:]
#    axFltArea.plot(yr, areaFlt, label=name)
#
#    GLflux = f.variables['groundingLineFlux'][:]
#    axGLflux.plot(yr, GLflux, label=name)
#
#    calvFlux = f.variables['totalCalvingFlux'][:]
#    axCalvFlux.plot(yr, calvFlux, label=name)

    # Fig 2: Grd MB ------------
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, volGround[:,r] - volGround[0,r], label="vol chg", linestyle=sty, color='k', linewidth=2)

    grdSMB = f.variables['regionalSumGroundedSfcMassBal'][:]
    if options.units == "m3":
       grdSMB = grdSMB / rhoi
    elif options.units == "kg":
       pass
    elif options.units == "Gt":
       grdSMB = grdSMB / 1.0e12
    cumGrdSMB = np.cumsum(grdSMB, axis=0)
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, cumGrdSMB[:,r], label="SMB", linestyle=sty, color='b')

    GLflux = f.variables['regionalSumGroundingLineFlux'][:]
    if options.units == "m3":
       GLflux = GLflux / rhoi
    elif options.units == "kg":
       pass
    elif options.units == "Gt":
       GLflux = GLflux / 1.0e12
    cumGLflux = np.cumsum(GLflux, axis=0)
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, -1.0*cumGLflux[:,r], label="GL flux", linestyle=sty, color='g')

    GLMigflux = f.variables['regionalSumGroundingLineMigrationFlux'][:]
    if options.units == "m3":
       GLMigflux = GLMigflux / rhoi
    elif options.units == "kg":
       pass
    elif options.units == "Gt":
       GLMigflux = GLMigflux / 1.0e12
    cumGLMigflux = np.cumsum(GLMigflux, axis=0)
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, -1.0*cumGLMigflux[:,r], label="GL mig flux", linestyle=sty, color='y')

    # sum of components
    grdSum = grdSMB - GLflux - GLMigflux # note negative sign on two GL terms - they are both positive grounded to floating
    cumGrdSum = np.cumsum(grdSum, axis=0)
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, cumGrdSum[:,r], label="sum", linestyle=sty, color='hotpink', linewidth=0.75)

    # Fig 3: Flt MB ---------------
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, volFloat[:,r] - volFloat[0,r], label="vol chg", linestyle=sty, color='k', linewidth=2)

    fltSMB = f.variables['regionalSumFloatingSfcMassBal'][:]
    if options.units == "m3":
       fltSMB = fltSMB / rhoi
    elif options.units == "kg":
       pass
    elif options.units == "Gt":
       fltSMB = fltSMB / 1.0e12
    cumFltSMB = np.cumsum(fltSMB, axis=0)
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumFltSMB[:,r], label="SMB", linestyle=sty, color='b')

    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumGLflux[:,r], label="GL flux", linestyle=sty, color='g')

    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumGLMigflux[:,r], label="GL mig flux", linestyle=sty, color='y')

    clv = f.variables['regionalSumCalvingFlux'][:]
    if options.units == "m3":
       clv = clv / rhoi
    elif options.units == "kg":
       pass
    elif options.units == "Gt":
       clv = clv / 1.0e12
    cumClv = np.cumsum(clv, axis=0)
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, -1.0*cumClv[:,r], label="calving", linestyle=sty, color='m', linewidth=1)

    BMB = f.variables['regionalSumFloatingBasalMassBal'][:]
    if options.units == "m3":
       BMB = BMB / rhoi
    elif options.units == "kg":
       pass
    elif options.units == "Gt":
       BMB = BMB / 1.0e12
    cumBMB = np.cumsum(BMB, axis=0)
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumBMB[:,r], label="BMB", linestyle=sty, color='r', linewidth=1)

    # sum of components
    fltSum = fltSMB + GLflux + GLMigflux - clv + BMB
    cumFltSum = np.cumsum(fltSum, axis=0)
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumFltSum[:,r], label="sum", linestyle=sty, color='hotpink', linewidth=0.75)




    f.close()


plotStat(options.file1inName, sty='-')


if(options.file2inName):
    plotStat(options.file2inName, sty=':')

if(options.file3inName):
    plotStat(options.file3inName, sty='--')

if(options.file4inName):
    plotStat(options.file4inName, sty='-.')


axs1.flatten()[-1].legend(loc='best', prop={'size': 6})
axs2.flatten()[-1].legend(loc='best', prop={'size': 6})
axs3.flatten()[-1].legend(loc='best', prop={'size': 6})

print("Generating plot.")
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
plt.show()

