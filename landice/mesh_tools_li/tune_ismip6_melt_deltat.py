#!/usr/bin/env python
'''
Script to calculate deltat for the ISMIP6 melt param to best match observed melt rates
for Antarctica.

Note that gamma0 is set in the script, as well as the range of deltaT to search over.

The tuned basin-by-basin deltaT values are written to a file called
'basin_and_coeff_gamma0_DeltaT_quadratic_non_local_gammaX.nc'. You will want to rename it to
avoid clobbering it if the script is rerun.  In addition to saving deltaT, the file also
includes the basin info and gamma0 so it can dropped directly into a streams file to run with.

For high res meshes, there may be efficiency gains that can be implemented, or the code may
need to move to Fortran.

Note: Currently hardcoded to use the first time level in the thermal forcing file.

Matt Hoffman, 9/8/2022
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt


# -----------------------
# --- Set gamma0 here ---
gamma0 = 14500.0 # MeanAnt
#gamma0 = 159000.0 # PIGL
# -----------------------
# Select range of deltaT values to search through.
# np.arange(-1.0, 1.5, 0.05) has been wide enough for MeanAnt with default gamma0 values,
# but range will be affected by gamma0 and a wider range may be necessary if any optimal deltaTs are outside this range.
# Confirm that the output deltaTs are more than increment away from boundary.
# Increments of 0.05 seems than fine enough to get a smooth function for interpolating, but use larger increments
# when testing for faster execution.
#dTs = np.arange(-1.5, 2.0, 0.25)  # MeanAnt - coarse spacing for rapid testing
dTs = np.arange(-1.0, 1.5, 0.05)  # MeanAnt - fine spacing for accurate calculation
#dTs = np.arange(-1.5, 0.0, 0.05)  # PIGL
# -----------------------


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-g", dest="fileName", help="input filename that includes ice geometry information.", metavar="FILENAME")
parser.add_option("-n", dest="fileRegionNames", help="region name filename.", metavar="FILENAME")
parser.add_option("-o", dest="fcgFileName", help="ocean forcing filename.  uses first time level", metavar="FILENAME")
options, args = parser.parse_args()

# Antarctic data from:
# Rignot, E., Bamber, J., van den Broeke, M. et al. Recent Antarctic ice mass loss from radar interferometry
# and regional climate modelling. Nature Geosci 1, 106-110 (2008). https://doi.org/10.1038/ngeo102
# Table 1: Mass balance of Antarctica in gigatonnes (10^12 kg) per year by sector for the year 2000
# https://www.nature.com/articles/ngeo102/tables/1
# and
# Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl. 2013. Ice-Shelf Melting Around Antarctica. Science 341 (6143): 266-70. https://doi.org/10.1126/science.1235798.
# Note: Only basin names and shelfMelt fields used in this script.
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

# Get region names from file
fn = Dataset(options.fileRegionNames, 'r')
rNamesIn = fn.variables['regionNames'][:]
regionCellMasks = fn.variables['regionCellMasks'][:]
nRegions = len(fn.dimensions['nRegions'])
# Process region names
rNamesOrig = list()
for reg in range(nRegions):
    thisString = rNamesIn[reg, :].tobytes().decode('utf-8').strip()  # convert from char array to string
    rNamesOrig.append(''.join(filter(str.isalnum, thisString)))  # this bit removes non-alphanumeric chars

# Parse region names to more usable names, if available
rNames = [None]*nRegions
for reg in range(nRegions):
    if rNamesOrig[reg] in ISMIP6basinInfo:
        rNames[reg] = ISMIP6basinInfo[rNamesOrig[reg]]['name']
    else:
        rNames[reg] = rNamesOrig[reg]

ff = Dataset(options.fcgFileName, 'r')
zOcean = ff.variables['ismip6shelfMelt_zOcean'][:]
TFocean = ff.variables['ismip6shelfMelt_3dThermalForcing'][0,:,:]

rhoi = 910.0
rhosw = 1028.0
cp_seawater = 3.974e3
latent_heat_ice = 335.0e3

coef = (rhosw*cp_seawater/(rhoi*latent_heat_ice))**2  # in K^(-2)

f = Dataset(options.fileName,'r')
thickness = f.variables['thickness'][0,:]
bedTopography = f.variables['bedTopography'][0,:]
floatMask = ((thickness*910/1028+bedTopography)<0)*(thickness>0)
lowerSurface = -rhoi/rhosw*thickness  #only works for floating areas
nCells = len(f.dimensions['nCells'])
areaCell = f.variables['areaCell'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]

def calcMelt(deltaT):
    TFdraft = np.zeros((nCells,))
    meanTFcell = np.zeros((nCells,))

    for iCell in range(nCells):
        if floatMask[iCell] == 1:
           # Linear interpolation of the thermal forcing on the ice draft depth:
           TFdraft[iCell] = np.interp(lowerSurface[iCell], np.flip(zOcean), np.flip(TFocean[iCell,:])) # flip b/c z is ordered from sfc to deep but interp needs to be increasing
    meanTF = np.zeros((nRegions,))
    for reg in range(nRegions):
        ind = np.nonzero(floatMask * regionCellMasks[:,reg])[0]
        meanTF[reg] = (TFdraft[ind]*areaCell[ind]).sum() / areaCell[ind].sum()
        meanTFcell[ind] = meanTF[reg]

    melt = gamma0 * coef * (TFdraft + deltaT) * np.absolute(meanTFcell + deltaT) # m/yr

    totalMelt = np.zeros((nRegions,))
    for reg in range(nRegions):
        ind = np.nonzero(floatMask * regionCellMasks[:,reg])[0]
        totalMelt[reg] = (melt[ind]*areaCell[ind]).sum() * rhoi / 1.0e12 #convert to Gt/yr

    return melt, totalMelt, TFdraft, meanTF

print("Considering dTs", dTs)
allMelts = np.zeros((len(dTs), nRegions))
for i, dT in enumerate(dTs):
    print(f"Calculating with dT={dT}")
    melt, allMelts[i,:], TFdraft, meanTF = calcMelt(dT)
    if False: # Generally don't want to produce a melt map for every dT, but this can be useful for debugging / special cases
       fig, axs = plt.subplots(1,1, figsize=(9,9), num=i)
       plt.scatter(xCell, yCell, s=1, c=melt, cmap='jet')
       plt.title(f'dT={dT}')
       plt.colorbar()

# Find dT that results in Rignot et al. (2013) melt rate for each basin
bestdTCells = np.zeros((nCells,))
regionCells = np.zeros((nCells,), 'i')
bestdT = np.zeros((nRegions,))
for reg in range(nRegions):
    bestdT[reg] = np.interp(ISMIP6basinInfo[rNamesOrig[reg]]['shelfMelt'][0], allMelts[:,reg], dTs)
    print(f"{ISMIP6basinInfo[rNamesOrig[reg]]['name']}: {bestdT[reg]}")
    bestdTCells[regionCellMasks[:,reg]==1] = bestdT[reg]
    # Also write out a region mask.
    # Note that regionCellMasks has a separate 0/1 mask for each region, whereas the ISMIP6 region mask is a single integer field where each cell is marked with the region to which it belongs.
    regionCells[regionCellMasks[:,reg]==1] = reg+1
#np.save(f'standard_bestdTs_{gamma0}.npy', bestdT)

nrow=4
ncol=4
fig4, axs4 = plt.subplots(nrow, ncol, figsize=(13, 11), num=4)
fig4.suptitle(f'melt sensitivity, gamma={gamma0}')
for reg in range(nRegions):
   plt.sca(axs4.flatten()[reg])
   plt.xlabel('delta T')
   plt.ylabel('total basin melt (Gt)')
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs4.flatten()[reg].set_title(f'{rNames[reg]}: {bestdT[reg]:.5f}')
   if reg == 0:
      axX = axs4.flatten()[reg]
   else:
      axs4.flatten()[reg].sharex(axX)
   meltHere = ISMIP6basinInfo[rNamesOrig[reg]]['shelfMelt'][0]
   plt.plot(dTs, meltHere * np.ones(dTs.shape), 'k-')
   plt.plot(dTs, allMelts[:,reg], 'b-')
fig4.tight_layout()

fig5, axs5 = plt.subplots(nrow, ncol, figsize=(13, 11), num=5)
fig5.suptitle(f'melt sensitivity, gamma={gamma0}')
for reg in range(nRegions):
   plt.sca(axs5.flatten()[reg])
   plt.xlabel('mean TF')
   plt.ylabel('total basin melt (Gt)')
   #plt.xticks(np.arange(22)*xtickSpacing)
   plt.grid()
   axs5.flatten()[reg].set_title(f'{rNames[reg]}: best TF={bestdT[reg]+meanTF[reg]:.3f}')
   if reg == 0:
      axX = axs5.flatten()[reg]
   else:
      axs5.flatten()[reg].sharex(axX)
   meltHere = ISMIP6basinInfo[rNamesOrig[reg]]['shelfMelt'][0]
   plt.plot(dTs+meanTF[reg], meltHere * np.ones(dTs.shape), 'k-')
   plt.plot(dTs+meanTF[reg] - bestdT[reg], allMelts[:,reg], 'b-')
   plt.plot(meanTF[reg]*np.ones((2,)), [0,meltHere*2.0], 'r--')
   
   TFs = dTs+meanTF[reg] #+ bestdTstd[reg]
   ind = np.nonzero(floatMask * regionCellMasks[:,reg])[0]
   plt.plot(TFs, coef * gamma0 * (TFs**2 + 2.0*bestdT[reg]*TFs + bestdT[reg]**2) * areaCell[ind].sum() * rhoi / 1.0e12, 'm:')
   
   #np.save(f'standard_allMelts_{gamma0}_{reg}.npy', allMelts)
fig5.tight_layout()
#np.save(f'meanTF.npy', meanTF)


# write new file
foutName=f'basin_and_coeff_DeltaT_quadratic_non_local_gamma{int(gamma0)}.nc'
fout = Dataset(foutName, 'w')
fout.createDimension('nCells', nCells)
dTOut = fout.createVariable('ismip6shelfMelt_deltaT', 'd', ('nCells',))
basinOut = fout.createVariable('ismip6shelfMelt_basin', 'i', ('nCells',))
gammaOut = fout.createVariable('ismip6shelfMelt_gamma0', 'd')
dTOut[:] = bestdTCells
basinOut[:] = regionCells
gammaOut[:] = gamma0
fout.close()

plt.show()
