#!/usr/bin/env python
'''
Script to plot calving information for a small ensemble of runs used for calving calibration.
The typical application is to run a handful of identical runs with different values of the
von Mises threshold stress.  The resulting plots can be used to assess which value gives
results that are stable and/or most similar to a control run with a fixed calving front
position, on a region by region basis.

It will look for files names regionalStats.nc in all subdirectories of the directory from
which it is run.  A control run with a fixed calving front should also be specified for
reference.

Matt Hoffman, Fall 2022
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import glob

rhoi = 910.0


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="Gt", metavar="FILENAME")
parser.add_option("-n", dest="fileRegionNames", help="region name filename.  If not specified, will attempt to read region names from file 1.", metavar="FILENAME")
parser.add_option("-c", dest="fileControl", help="control run with fixed calving front", metavar="FILENAME")
options, args = parser.parse_args()

print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))

runs = sorted(glob.glob('*/regionalStats.nc'))
print(runs)
nRuns = len(runs)

if options.units == "m3":
   massUnit = "m$^3$"
elif options.units == "kg":
   massUnit = "kg"
elif options.units == "Gt":
   massUnit = "Gt"
else:
   sys.exit("Unknown mass/volume units")
print("Using volume/mass units of: ", massUnit)

# Get region names from file
fn = Dataset(options.fileRegionNames, 'r')
nRegions = len(fn.dimensions['nRegions'])
rNamesIn = fn.variables['regionNames'][:]

# Process region names
rNamesOrig = list()
for r in range(nRegions):
    thisString = rNamesIn[r, :].tobytes().decode('utf-8').strip()  # convert from char array to string
    rNamesOrig.append(''.join(filter(str.isalnum, thisString)))  # this bit removes non-alphanumeric chars

# Antarctic data from:
# Rignot, E., Bamber, J., van den Broeke, M. et al. Recent Antarctic ice mass loss from radar interferometry
# and regional climate modelling. Nature Geosci 1, 106-110 (2008). https://doi.org/10.1038/ngeo102
# Table 1: Mass balance of Antarctica in gigatonnes (10^12 kg) per year by sector for the year 2000
# https://www.nature.com/articles/ngeo102/tables/1
# Note: May want to switch to input+, net+
# Note: Some ISMIP6 basins combine multiple Rignot basins.  May want to separate if we update our regions.
ISMIP6basinInfo = {
        'ISMIP6BasinAAp': {'name': 'Dronning Maud Land', 'input': [60,9], 'outflow': [60,7], 'net': [0, 11], 'shelfMelt': [57.5]},
        'ISMIP6BasinApB': {'name': 'Enderby Land', 'input': [39,5], 'outflow': [40,2], 'net': [-1,5], 'shelfMelt': [24.6]},
        'ISMIP6BasinBC': {'name': 'Amery-Lambert', 'input': [73, 10], 'outflow': [77,4], 'net': [-4, 11], 'shelfMelt': [35.5]},
        'ISMIP6BasinCCp': {'name': 'Phillipi, Denman', 'input': [81, 13], 'outflow': [87,7], 'net':[-7,15], 'shelfMelt': [107.9]},
        'ISMIP6BasinCpD': {'name': 'Totten', 'input': [198,37], 'outflow': [207,13], 'net': [-8,39], 'shelfMelt': [102.3]},
        'ISMIP6BasinDDp': {'name': 'Mertz', 'input': [93,14], 'outflow': [94,6], 'net': [-2,16], 'shelfMelt': [22.8]},
        'ISMIP6BasinDpE': {'name': 'Victoria Land', 'input': [20,1], 'outflow': [22,3], 'net': [-2,4], 'shelfMelt': [22.9]},
        'ISMIP6BasinEF': {'name': 'Ross', 'input': [61+110,(10**2+7**2)**0.5], 'outflow': [49+80,(4**2+2^2)**0.5], 'net': [11+31,(11*2+7**2)**0.5], 'shelfMelt': [70.3]},
        'ISMIP6BasinFG': {'name': 'Getz', 'input': [108,28], 'outflow': [128,18], 'net': [-19,33], 'shelfMelt': [152.9]},
        'ISMIP6BasinGH': {'name': 'Thwaites/PIG', 'input': [177,25], 'outflow': [237,4], 'net': [-61,26], 'shelfMelt': [290.9]},
        'ISMIP6BasinHHp': {'name': 'Bellingshausen', 'input': [51,16], 'outflow': [86,10], 'net': [-35,19], 'shelfMelt': [76.3]},
        'ISMIP6BasinHpI': {'name': 'George VI', 'input': [71,21], 'outflow': [78,7], 'net': [-7,23], 'shelfMelt': [152.3]},
        'ISMIP6BasinIIpp': {'name': 'Larsen A-C', 'input': [15,5], 'outflow': [20,3], 'net': [-5,6], 'shelfMelt': [32.9]},
        'ISMIP6BasinIppJ': {'name': 'Larsen E', 'input': [8,4], 'outflow': [9,2], 'net': [-1,4], 'shelfMelt': [4.3]},
        'ISMIP6BasinJK': {'name': 'FRIS', 'input': [93+142, (8**2+11**2)**0.5], 'outflow': [75+145,(4**2+7**2)**0.5], 'net': [18-4,(9**2+13**2)**0.5], 'shelfMelt': [155.4]},
        'ISMIP6BasinKA': {'name': 'Brunt-Stancomb', 'input': [42+26,(8**2+7**2)**0.5], 'outflow': [45+28,(4**2+2**2)**0.5], 'net':[-3-1,(9**2+8**2)**0.5], 'shelfMelt': [10.4]}
        }

# Parse region names to more usable names, if available
rNames = [None]*nRegions
for r in range(nRegions):
    if rNamesOrig[r] in ISMIP6basinInfo:
        rNames[r] = ISMIP6basinInfo[rNamesOrig[r]]['name']
    else:
        rNames[r] = rNamesOrig[r]

#print(rNames)

nrow=4
ncol=4
if nRegions > nrow*ncol:
    sys.exit("ERROR: Number of regions exceeds number of plots.  Please adjust nrow, ncol as needed.")

# Set up Figure 1: total area
fig1, axs1 = plt.subplots(nrow, ncol, figsize=(13, 11), num=1)
fig1.suptitle(f'area change')
for reg in range(nRegions):
   plt.sca(axs1.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('area change (km^2)')
   plt.grid()
   axs1.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs1.flatten()[reg]
   else:
      axs1.flatten()[reg].sharex(axX)

# Set up Figure 2: grd vol
fig2, axs2 = plt.subplots(nrow, ncol, figsize=(13, 11), num=2)
fig2.suptitle(f'grounded volume change')
for reg in range(nRegions):
   plt.sca(axs2.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('volume change ({})'.format(massUnit))
   plt.grid()
   axs2.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs2.flatten()[reg]
   else:
      axs2.flatten()[reg].sharex(axX)

# Set up Figure 3: VAF
fig3, axs3 = plt.subplots(nrow, ncol, figsize=(13, 11), num=3)
fig3.suptitle(f'VAF change')
for reg in range(nRegions):
   plt.sca(axs3.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('VAF change ({})'.format(massUnit))
   plt.grid()
   axs3.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs3.flatten()[reg]
   else:
      axs3.flatten()[reg].sharex(axX)

# Set up Figure 4: calving
fig4, axs4 = plt.subplots(nrow, ncol, figsize=(13, 11), num=4)
fig4.suptitle(f'Calving rate')
for reg in range(nRegions):
   plt.sca(axs4.flatten()[reg])
   if reg // nrow == nrow-1:
      plt.xlabel('Year')
   if reg % ncol == 0:
      plt.ylabel('Calving rate ({}/yr)'.format(massUnit))
   plt.grid()
   axs4.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs4.flatten()[reg]
   else:
      axs4.flatten()[reg].sharex(axX)


def plotStat(fname, col, addToLegend=False):
    print("Reading and plotting file: {}".format(fname))

    sty='-'

    name = fname
    runName = fname.split('/')[0]

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    dt = f.variables['deltat'][:]/(3600.0*24.0*365.0) # in yr
    dtnR = np.tile(dt.reshape(len(dt),1), (1,nRegions))  # repeated per region with dim of nt,nRegions
    nRegionsLocal = len(f.dimensions['nRegions'])
    if nRegionsLocal != nRegions:
        sys.exit(f"ERROR: Number of regions in file {fname} does not match number of regions in first input file!")

    # Fig 1: area change  ---------------
    areaTot = f.variables['regionalIceArea'][:]/1000.0**2
    areaGrd = f.variables['regionalGroundedIceArea'][:]/1000.0**2
    areaFlt = f.variables['regionalFloatingIceArea'][:]/1000.0**2
    for r in range(nRegions):
        axs1.flatten()[r].plot(yr, areaTot[:,r] - areaTot[0,r], label=f"{runName}", linestyle=sty, color=col)

    # Fig 2: grd vol change  ---------------
    volGround = f.variables['regionalGroundedIceVolume'][:]
    volGround = volGround[:,:] - volGround[0,:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       volGround = volGround * rhoi
    elif options.units == "Gt":
       volGround = volGround * rhoi / 1.0e12
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, volGround[:,r] - volGround[0,r], label=f"{runName}", linestyle=sty, color=col)

    # Fig 3: VAF change  ---------------
    VAF = f.variables['regionalVolumeAboveFloatation'][:]
    VAF = VAF[:,:] - VAF[0,:]
    if options.units == "m3":
       pass
    elif options.units == "kg":
       VAF = VAF * rhoi
    elif options.units == "Gt":
       VAF = VAF * rhoi / 1.0e12
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, VAF[:,r] - VAF[0,r], label=f"{runName}", linestyle=sty, color=col)

    # Fig 4: calving ------------
    clv = f.variables['regionalSumCalvingFlux'][:]
    if options.units == "m3":
       clv = clv / rhoi
    elif options.units == "kg":
       pass
    elif options.units == "Gt":
       clv = clv / 1.0e12
    cumClv = np.cumsum(clv*dtnR, axis=0)
    lbl = 'calving' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs4.flatten()[r].plot(yr, clv[:,r], label=f"{runName}", linestyle=sty, color=col)


    f.close()


plotStat(options.fileControl, col='k', addToLegend=True)
colors = pl.cm.jet(np.linspace(0,1,nRuns))
for cnt, run in enumerate(runs):
    plotStat(run, col=colors[cnt], addToLegend=True)


axs1.flatten()[-1].legend(loc='best', prop={'size': 5})
axs2.flatten()[-1].legend(loc='best', prop={'size': 5})
axs3.flatten()[-1].legend(loc='best', prop={'size': 5})
axs4.flatten()[-1].legend(loc='best', prop={'size': 5})

print("Generating plot.")
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
plt.show()

