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
parser.add_option("-n", dest="fileRegionNames", help="region name filename.  If not specified, will attempt to read region names from file 1.", metavar="FILENAME")
options, args = parser.parse_args()

print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))

# Build string for titles about the runs in use
runinfo=f'solid={options.file1inName}'
if options.file2inName:
    runinfo = f'{runinfo}\ndotted={options.file2inName}'
if options.file3inName:
    runinfo = f'{runinfo}\ndashed={options.file3inName}'
if options.file4inName:
    runinfo = f'{runinfo}\ndashdot={options.file4inName}'

#xtickSpacing = 20.0

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

# Get nRegions and yr from first file
f = Dataset(options.file1inName, 'r')
nRegions = len(f.dimensions['nRegions'])
yr = f.variables['daysSinceStart'][:]/365.0

# Get region names from file
if options.fileRegionNames:
   fn = Dataset(options.fileRegionNames, 'r')
   rNamesIn = fn.variables['regionNames'][:]
else:
   rNamesIn = f.variables['regionNames'][:]
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
# and
# Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl. 2013. Ice-Shelf Melting Around Antarctica. Science 341 (6143): 266-70. https://doi.org/10.1126/science.1235798.
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

if nRegions <= 4:
    ncol = 2
elif nRegions <= 9:
    ncol = 3
elif nRegions <= 16:
    ncol = 4
elif nRegions <= 25:
    ncol = 5
else:
    sys.exit("ERROR: More than 25 regions found.  Attempting to plot this many regions is likely a bad idea.")
nrow = np.ceil(nRegions / ncol).astype('int') # Set nrow to have enough rows to plot number of regions based on ncol calculated above

# Set up Figure 1: volume stats overview
fig1, axs1 = plt.subplots(nrow, ncol, figsize=(13, 11), num=1)
fig1.suptitle(f'Mass change summary\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs1.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('volume change ({})'.format(massUnit))
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs1.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs1.flatten()[reg]
   else:
      axs1.flatten()[reg].sharex(axX)
   # plot obs if applicable
   if rNamesOrig[reg] in ISMIP6basinInfo:
       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['net']
       axs1.flatten()[reg].fill_between(yr, yr*(mn-sig), yr*(mn+sig), color='b', alpha=0.2, label='grd obs')

# Set up Figure 2: grounded MB
fig2, axs2 = plt.subplots(nrow, ncol, figsize=(13, 11), num=2)
fig2.suptitle(f'Grounded mass change\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs2.flatten()[reg])
   if reg // nrow == nrow-1:
      plt.xlabel('Year')
   if reg % ncol == 0:
      plt.ylabel('volume change ({})'.format(massUnit))
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs2.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs2.flatten()[reg]
   else:
      axs2.flatten()[reg].sharex(axX)
   # plot obs if applicable
   if rNamesOrig[reg] in ISMIP6basinInfo:
       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['input']
       axs2.flatten()[reg].fill_between(yr, yr*(mn-sig), yr*(mn+sig), color='b', alpha=0.2, label='SMB obs')
       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['outflow']
       axs2.flatten()[reg].fill_between(yr, -yr*(mn-sig), -yr*(mn+sig), color='g', alpha=0.2, label='outflow obs')
       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['net']
       axs2.flatten()[reg].fill_between(yr, yr*(mn-sig), yr*(mn+sig), color='k', alpha=0.2, label='net obs')


# Set up Figure 3: floating MB
fig3, axs3 = plt.subplots(nrow, ncol, figsize=(13, 11), num=3)
fig3.suptitle(f'Floating mass change\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs3.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('volume change ({})'.format(massUnit))
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs3.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs3.flatten()[reg]
   else:
      axs3.flatten()[reg].sharex(axX)

# Set up Figure 4: area change
fig4, axs4 = plt.subplots(nrow, ncol, figsize=(13, 11), num=4)
fig4.suptitle(f'Area change\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs4.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('Area change (km^2)')
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs4.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs4.flatten()[reg]
   else:
      axs4.flatten()[reg].sharex(axX)


# Set up Figure 5
fig5, axs5 = plt.subplots(2,1, figsize=(13, 11), num=5)
fig5.suptitle(f'regional contributions\n{runinfo}', fontsize=9)
mnTot=0.0
sigTot = 0.0
for reg in range(nRegions):
   if rNamesOrig[reg] in ISMIP6basinInfo:
       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['net']
       mnTot += mn
       sigTot += sig**2
sigTot = sigTot**0.5
axs5.flatten()[0].fill_between(yr, yr*(mnTot-sigTot), yr*(mnTot+sigTot), color='k', alpha=0.2, label='net obs')
plt.sca(axs5.flatten()[0])
plt.xlabel('Year')
plt.ylabel('Mass change (Gt)')
plt.grid()
axs5.flatten()[1].fill_between(yr, yr*(mnTot-sigTot), yr*(mnTot+sigTot), color='k', alpha=0.2, label='net obs')
plt.sca(axs5.flatten()[1])
plt.xlabel('Year')
plt.ylabel('VAF mass change (Gt)')
plt.grid()


# Set up Figure 6: melt rate vs obs
fig6, axs6 = plt.subplots(nrow, ncol, figsize=(13, 11), num=6)
fig6.suptitle(f'Ice-shelf melt rate\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs6.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('Ice-shelf melt rate (Gt/yr)')
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs6.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs6.flatten()[reg]
   else:
      axs6.flatten()[reg].sharex(axX)
   if rNamesOrig[reg] in ISMIP6basinInfo:
       mlt = ISMIP6basinInfo[rNamesOrig[reg]]['shelfMelt'][0]
       axs6.flatten()[reg].plot(yr, np.ones(yr.shape)*(mlt), color='k', label='melt obs')

# Set up unit conversion factors to be used when reading variables
if options.units == "m3":
    volUnitFactor = 1.0
    massUnitFactor = 1.0 / rhoi
elif options.units == "kg":
    volUnitFactor = rhoi
    massUnitFactor = 1.0
elif options.units == "Gt":
    volUnitFactor = rhoi / 1.0e12
    massUnitFactor = 1.0 / 1.0e12
else:
    sys.exit("ERROR: Unknown unit specified")


def plotStat(fname, sty, addToLegend=False):
    print("Reading and plotting file: {}".format(fname))

    name = fname

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    dt = f.variables['deltat'][:]/(3600.0*24.0*365.0) # in yr
    #yr = yr-yr[0]  # uncomment to align all start dates
    dtnR = np.tile(dt.reshape(len(dt),1), (1,nRegions))  # repeated per region with dim of nt,nRegions
    nRegionsLocal = len(f.dimensions['nRegions'])
    if nRegionsLocal != nRegions:
        sys.exit(f"ERROR: Number of regions in file {fname} does not match number of regions in first input file!")

    # Fig 1: summary plot
    vol = f.variables['regionalIceVolume'][:] * volUnitFactor
    lbl ='total' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, vol[:,r] - vol[0,r], label=lbl, linestyle=sty, color='k')

    VAF = f.variables['regionalVolumeAboveFloatation'][:] * volUnitFactor
    VAF = VAF[:,:] - VAF[0,:]
    lbl ='VAF' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, VAF[:,r] - VAF[0,r], label=lbl, linestyle=sty, color='m')

    volGround = f.variables['regionalGroundedIceVolume'][:] * volUnitFactor
    volGround = volGround[:,:] - volGround[0,:]
    lbl ='grd' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, volGround[:,r] - volGround[0,r], label=lbl, linestyle=sty, color='b')

    volFloat = f.variables['regionalFloatingIceVolume'][:] * volUnitFactor
    volFloat = volFloat[:,:] - volFloat[0,:]
    lbl ='flt' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, volFloat[:,r] - volFloat[0,r], label=lbl, linestyle=sty, color='g')


    # Fig 2: Grd MB ------------
    lbl ='vol chg' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, volGround[:,r] - volGround[0,r], label=lbl, linestyle=sty, color='k', linewidth=2)

    grdSMB = f.variables['regionalSumGroundedSfcMassBal'][:] * massUnitFactor
    cumGrdSMB = np.cumsum(grdSMB*dtnR, axis=0)
    lbl ='SMB' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, cumGrdSMB[:,r], label=lbl, linestyle=sty, color='b')

    GLflux = f.variables['regionalSumGroundingLineFlux'][:] * massUnitFactor
    cumGLflux = np.cumsum(GLflux*dtnR, axis=0)
    lbl ='GL flux' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, -1.0*cumGLflux[:,r], label=lbl, linestyle=sty, color='g')

    GLMigflux = f.variables['regionalSumGroundingLineMigrationFlux'][:] * massUnitFactor
    cumGLMigflux = np.cumsum(GLMigflux*dtnR, axis=0)
    lbl ='GL mig flux' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, -1.0*cumGLMigflux[:,r], label=lbl, linestyle=sty, color='y')

    # sum of components
    grdSum = grdSMB - GLflux - GLMigflux # note negative sign on two GL terms - they are both positive grounded to floating
    cumGrdSum = np.cumsum(grdSum*dtnR, axis=0)
    lbl ='sum' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, cumGrdSum[:,r], label=lbl, linestyle=sty, color='hotpink', linewidth=0.75)
    grdSum2 = grdSMB - GLflux  # note negative sign on two GL terms - they are both positive grounded to floating
    cumGrdSum2 = np.cumsum(grdSum2*dtnR, axis=0)
    lbl ='sum, no GLmig' if addToLegend else '_nolegend_'
    for r in range(nRegions):
        axs2.flatten()[r].plot(yr, cumGrdSum2[:,r], label=lbl, linestyle=':', color='hotpink', linewidth=0.75)


    # Fig 3: Flt MB ---------------
    lbl ='vol chg' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, volFloat[:,r] - volFloat[0,r], label=lbl, linestyle=sty, color='k', linewidth=2)

    fltSMB = f.variables['regionalSumFloatingSfcMassBal'][:] * massUnitFactor
    cumFltSMB = np.cumsum(fltSMB*dtnR, axis=0)
    lbl = 'SMB' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumFltSMB[:,r], label=lbl, linestyle=sty, color='b')

    lbl = 'GL flux' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumGLflux[:,r], label=lbl, linestyle=sty, color='g')

    lbl = 'GL mig flux' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumGLMigflux[:,r], label=lbl, linestyle=sty, color='y')

    clv = f.variables['regionalSumCalvingFlux'][:] * massUnitFactor
    cumClv = np.cumsum(clv*dtnR, axis=0)
    lbl = 'calving' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, -1.0*cumClv[:,r], label=lbl, linestyle=sty, color='m', linewidth=1)

    BMB = f.variables['regionalSumFloatingBasalMassBal'][:] * massUnitFactor
    cumBMB = np.cumsum(BMB*dtnR, axis=0)
    lbl = 'BMB' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumBMB[:,r], label=lbl, linestyle=sty, color='r', linewidth=1)

    # sum of components
    fltSum = fltSMB + GLflux + GLMigflux - clv + BMB
    cumFltSum = np.cumsum(fltSum*dtnR, axis=0)
    lbl = 'sum' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumFltSum[:,r], label=lbl, linestyle=sty, color='hotpink', linewidth=0.75)
    fltSum2 = fltSMB + GLflux - clv + BMB
    cumFltSum2 = np.cumsum(fltSum2*dtnR, axis=0)
    lbl = 'sum, no GLmig' if addToLegend else '_nolegend_'
    for r in range(nRegions):
        axs3.flatten()[r].plot(yr, cumFltSum2[:,r], label=lbl, linestyle=':', color='hotpink', linewidth=0.75)


    # Fig 4: area change  ---------------

    areaTot = f.variables['regionalIceArea'][:]/1000.0**2
    areaGrd = f.variables['regionalGroundedIceArea'][:]/1000.0**2
    areaFlt = f.variables['regionalFloatingIceArea'][:]/1000.0**2
    for r in range(nRegions):
        axs4.flatten()[r].plot(yr, areaTot[:,r] - areaTot[0,r], label=("total area" if addToLegend else '_nolegend_'), linestyle=sty, color='k')
        axs4.flatten()[r].plot(yr, areaGrd[:,r] - areaGrd[0,r], label=("grd area" if addToLegend else '_nolegend_'), linestyle=sty, color='b')
        axs4.flatten()[r].plot(yr, areaFlt[:,r] - areaFlt[0,r], label=("flt area" if addToLegend else '_nolegend_'), linestyle=sty, color='g')

    # Fig. 5:  select global stats ---------
    for r in range(nRegions):
        if rNamesOrig[r] == 'ISMIP6BasinGH':
           indTG = r
           break
    #print(f'TG index={indTG}')
    axs5.flatten()[0].plot(yr, volGround.sum(axis=1), label='total', color='b', linestyle=sty)
    volGroundnoTG = np.delete(volGround, indTG, 1)
    axs5.flatten()[0].plot(yr, volGroundnoTG.sum(axis=1), label='no TG/PIG', color='c', linestyle=sty)
    #for r in range(nRegions):
       #axs5.flatten()[0].plot(yr, VAF[:,r] - VAF[0,r], label=rNames[r], linestyle=sty)
    axs5.flatten()[1].plot(yr, VAF.sum(axis=1), label='total', color='b', linestyle=sty)
    VAFnoTG = np.delete(VAF, indTG, 1)
    axs5.flatten()[1].plot(yr, VAFnoTG.sum(axis=1), label='no TG/PIG', color='c', linestyle=sty)

    # Fig. 6:  melt rates ---------
    #BMB = np.where(BMB==0.0, BMB, np.nan*BMB)
    for r in range(nRegions):
        axs6.flatten()[r].plot(yr, -BMB[:,r], label=("BMB" if addToLegend else '_nolegend_'), linestyle=sty, color='b')

    f.close()


plotStat(options.file1inName, sty='-', addToLegend=True)


if(options.file2inName):
    plotStat(options.file2inName, sty=':')

if(options.file3inName):
    plotStat(options.file3inName, sty='--')

if(options.file4inName):
    plotStat(options.file4inName, sty='-.')


axs1.flatten()[-1].legend(loc='best', prop={'size': 5})
axs2.flatten()[-1].legend(loc='best', prop={'size': 5})
axs3.flatten()[-1].legend(loc='best', prop={'size': 5})
axs4.flatten()[-1].legend(loc='best', prop={'size': 6})
axs5.flatten()[0].legend(loc='best', prop={'size': 6})
axs6.flatten()[-1].legend(loc='best', prop={'size': 6})

print("Generating plot.")
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()
fig6.tight_layout()
plt.show()

