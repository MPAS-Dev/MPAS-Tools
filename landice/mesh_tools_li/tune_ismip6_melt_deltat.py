#!/usr/bin/env python
'''
Script to calculate deltat for the ISMIP6 melt param to best match observed melt rates.
If 'shelfBaseSlope' is included in input file, the calculation will use the form of the
parameterization that includes slope.  To get this field, run MALI one timestep with
config_ismip6shelfMelt_use_slope=.true. and shelfBaseSlope in the output stream.  This will
ensure that the slope is calculated in an identical way to how MALI would do it.

Note that gamma0 is set in the script, as well as the range of deltaT to search over.

The tuned basin-by-basin deltaT values are written to a file called
'basin_and_coeff_gamma0_DeltaT_quadratic_non_local.nc'. You will want to rename it to
avoid clobbering it if the script is rerun.  In addition to saving deltaT, the file also
include the basin info and gamma0 so it can dropped directly into a streams file to run with.

For high res meshes, there may be efficiency gains that can be implemented, or the code may
need to move to Fortran.

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
#gamma0 = 14500.0 # MeanAnt
#gamma0 = 159000.0 # PIGL
gamma0 = 2.06e6 #MeanAnt Median with slope
# -----------------------
# Select range of deltaT values to search through.
# np.arange(-1.0, 1.5, 0.05) has been wide enough for MeanAnt with and without slope with default gamma0 values,
# but range will be affected by gamma0 and a wider range may be necessary if any optimal deltaTs are outside this range.
# Confirm that the output deltaTs are more than increment away from boundary.
# Increments of 0.05 seems than fine enough to get a smooth function for interpolating, but use larger increments
# when testing for faster execution.
dTs = np.arange(-1.0, 1.5, 0.05)
# -----------------------


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-g", dest="fileName", help="input filename that includes ice geometry information. If shelfBaseSlope is included in this file, the script will assume the form of the melt parameterization that includes slope", metavar="FILENAME")
parser.add_option("-n", dest="fileRegionNames", help="region name filename.", metavar="FILENAME")
parser.add_option("-o", dest="fcgFileName", help="ocean forcing filename", metavar="FILENAME")
options, args = parser.parse_args()

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
for r in range(nRegions):
    thisString = rNamesIn[r, :].tobytes().decode('utf-8').strip()  # convert from char array to string
    rNamesOrig.append(''.join(filter(str.isalnum, thisString)))  # this bit removes non-alphanumeric chars

# Parse region names to more usable names, if available
rNames = [None]*nRegions
for r in range(nRegions):
    if rNamesOrig[r] in ISMIP6basinInfo:
        rNames[r] = ISMIP6basinInfo[rNamesOrig[r]]['name']
    else:
        rNames[r] = rNamesOrig[r]

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
if 'shelfBaseSlope' in f.variables:
    print('USING SLOPE')
    slope = f.variables['shelfBaseSlope'][0,:]
else:
    slope = np.ones((nCells,))

def calcMelt(deltaT):
    mean_TF = 0.0
    IS_area = 0.0
    meanSlope = 0.0
    TFdraft = np.zeros((nCells,))
    meanTFcell = np.zeros((nCells,))

    for iCell in range(nCells):
        if floatMask[iCell] == 1:
           # 1 -  Linear interpolation of the thermal forcing on the ice draft depth :
           TFdraft[iCell] = np.interp(lowerSurface[iCell], np.flip(zOcean), np.flip(TFocean[iCell,:])) # flip b/c z is ordered from sfc to deep but interp needs to be increasing
    meanTF = np.zeros((nRegions,))
    for r in range(nRegions):
        ind = np.nonzero(floatMask * regionCellMasks[:,r])[0]
        meanTF[r] = (TFdraft[ind]*areaCell[ind]).sum() / areaCell[ind].sum()
        meanTFcell[ind] = meanTF[r]

    melt = gamma0 * coef * slope * (TFdraft + deltaT) * np.absolute(meanTFcell + deltaT) # m/yr

    totalMelt = np.zeros((nRegions,))
    for r in range(nRegions):
        ind = np.nonzero(floatMask * regionCellMasks[:,r])[0]
        totalMelt[r] = (melt[ind]*areaCell[ind]).sum() * rhoi / 1.0e12 #convert to Gt/yr

    return melt, totalMelt, TFdraft

print("Considering dTs", dTs)
allMelts = np.zeros((len(dTs), nRegions))
for i, dT in enumerate(dTs):
    print(f"Calculating with dT={dT}")
    melt, allMelts[i,:], TFdraft = calcMelt(dT)
    if False:
       fig, axs = plt.subplots(1,1, figsize=(9,9), num=i)
       plt.scatter(xCell, yCell, s=1, c=melt, cmap='jet')
       plt.title(f'dT={dT}')
       plt.colorbar()


bestdTCells = np.zeros((nCells,))
regionCells = np.zeros((nCells,), 'i')
bestdT = np.zeros((nRegions,))
for r in range(nRegions):
    bestdT[r] = np.interp(ISMIP6basinInfo[rNamesOrig[r]]['shelfMelt'][0], allMelts[:,r], dTs)
    print(f"{ISMIP6basinInfo[rNamesOrig[r]]['name']}: {bestdT[r]}")
    bestdTCells[regionCellMasks[:,r]==1] = bestdT[r]
    regionCells[regionCellMasks[:,r]==1] = r+1

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

# write new file
foutName='basin_and_coeff_gamma0_DeltaT_quadratic_non_local.nc'
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
