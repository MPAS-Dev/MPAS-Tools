#!/usr/bin/env python
'''
Script to plot subglacial hydrology time-series from a landice globalStats file.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import argparse

rhow = 1000.0
secyr = 3600.0 * 24.0 * 365.0

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-f", dest="filename", help="input filename", default="globalStats.nc", metavar="FILENAME")
parser.add_argument("-u", dest="units", help="units for mass: kg, Gt", default="Gt", metavar="FILENAME")
args = parser.parse_args()

# Scaling assuming variables are in kg
if args.units == "kg":
   massUnit = "kg"
   fluxUnit = "kg yr$^{-1}$"
   unitScaling = 1.0
elif args.units == "Gt":
   massUnit = "Gt"
   fluxUnit = "Gt yr$^{-1}$"
   unitScaling = 1.0e-12
else:
   sys.exit("Unknown mass unit")

print("Using mass units of: ", massUnit)


dataset = nc.Dataset(args.filename)

# Read variables
# convert everything to kg and years before unit conversion
totalSubglacialWaterMass = \
    dataset.variables['totalSubglacialWaterVolume'][:] * rhow * unitScaling
melt = dataset.variables['totalBasalMeltInput'][:] * unitScaling * secyr
distFluxMarine = dataset.variables['totalDistWaterFluxMarineMargin'][:] * unitScaling * secyr
chnlFluxMarine = dataset.variables['totalChnlWaterFluxMarineMargin'][:] * unitScaling * secyr
distFluxLand = dataset.variables['totalDistWaterFluxTerrestrialMargin'][:] * unitScaling * secyr
chnlFluxLand = dataset.variables['totalChnlWaterFluxTerrestrialMargin'][:] * unitScaling * secyr
chnlMelt = dataset.variables['totalChannelMelt'][:] * unitScaling * secyr
flotFrac = dataset.variables['avgFlotationFraction'][:]
lakeArea = dataset.variables['totalSubglacialLakeArea'][:] / 1000.0**2  # km^2
lakeMass = dataset.variables['totalSubglacialLakeVolume'][:] * rhow * unitScaling

deltat = dataset.variables['deltat'][:] / secyr  # in years
yr = dataset.variables['daysSinceStart'][:] / 365.0

subglacialWaterMassRate = np.zeros((len(melt),))

for i in range(len(totalSubglacialWaterMass) - 1):
    subglacialWaterMassRate[i] = ((totalSubglacialWaterMass[i+1] -
                                   totalSubglacialWaterMass[i]) / deltat[i])

# plot mass balance time-series
fig, ax = plt.subplots(1, 1, layout='tight', figsize=(8,6))
# input
plt.plot(yr, melt, 'r:', label='basal melt')
plt.plot(yr, chnlMelt, 'r--', label='channel melt')
total_melt = melt + chnlMelt
plt.plot(yr, total_melt, 'r-', label='total melt')

# output
plt.plot(yr, distFluxMarine, 'b--', label='marine sheet outflux')
plt.plot(yr, distFluxLand, 'b:', label='land sheet outflux')
plt.plot(yr, chnlFluxMarine, 'c--', label='marine chnl outflux')
plt.plot(yr, chnlFluxLand, 'c:', label='land chnl outflux')
total_outflux = distFluxMarine + distFluxLand + chnlFluxMarine + chnlFluxLand
plt.plot(yr, total_outflux, 'b-', lw=2, label='total outflux')

plt.plot(yr, subglacialWaterMassRate, 'g-', label='dV/dt')

plt.plot(yr, total_melt - total_outflux, 'k', label='I-O')

plt.legend(loc='best')
plt.xlabel('Year')
plt.ylabel(f'Mass flux ({fluxUnit})')

# Plot other time-series of interest
fig, axes = plt.subplots(2, 2, sharex=True, layout='tight',
                         figsize=(10, 7))
axes = axes.flatten()

ax = 0
axes[ax].plot(yr, flotFrac)
axes[ax].set_ylabel('Flotation fraction')

ax += 1
axes[ax].plot(yr, totalSubglacialWaterMass)
axes[ax].set_ylabel(f'Water mass ({massUnit})')

ax += 1
axes[ax].plot(yr, lakeArea)
axes[ax].set_ylabel('Lake area (km$^2$)')

ax += 1
axes[ax].plot(yr, lakeMass)
axes[ax].set_ylabel(f'Lake mass ({massUnit})')

for ax in axes:
    ax.grid(True)
    ax.set_xlabel("Year")


plt.show()

