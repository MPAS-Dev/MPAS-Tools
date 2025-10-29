#!/usr/bin/env python
'''
Script to plot subglacial hydrology time-series from a landice globalStats file.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import argparse

rhoi = 910.0
rhosw = 1028.

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-f", dest="filename", help="input filename", default="globalStats.nc", metavar="FILENAME")
args = parser.parse_args()

dataset = nc.Dataset(args.filename)


# plot mass balance time-series
totalSubglacialWaterVolume = dataset.variables['totalSubglacialWaterVolume'][:]
melt = dataset.variables['totalBasalMeltInput'][:]
distFluxMarine = dataset.variables['totalDistWaterFluxMarineMargin'][:]
chnlFluxMarine = dataset.variables['totalChnlWaterFluxMarineMargin'][:]
distFluxLand = dataset.variables['totalDistWaterFluxTerrestrialMargin'][:]
chnlFluxLand = dataset.variables['totalChnlWaterFluxTerrestrialMargin'][:]


chnlMelt = dataset.variables['totalChannelMelt'][:]

deltat = dataset.variables['deltat'][:]
yr = dataset.variables['daysSinceStart'][:] / 365.0

subglacialWaterMassRate = np.zeros((len(melt),))
rhow = 1000.0

for i in range(len(totalSubglacialWaterVolume) - 1):
    subglacialWaterMassRate[i] = ((totalSubglacialWaterVolume[i+1] - totalSubglacialWaterVolume[i]) * rhow / deltat[i])

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
plt.ylabel('Mass flux (kg/s)')

# Plot other time-series of interest
flotFrac = dataset.variables['avgFlotationFraction'][:]
lakeArea = dataset.variables['totalSubglacialLakeArea'][:]
lakeVol = dataset.variables['totalSubglacialLakeVolume'][:]

fig, axes = plt.subplots(2, 2, sharex=True, layout='tight',
                         figsize=(10, 7))
axes = axes.flatten()

ax = 0
axes[ax].plot(yr, flotFrac)
axes[ax].set_ylabel('Flotation fraction')

ax += 1
axes[ax].plot(yr, totalSubglacialWaterVolume)
axes[ax].set_ylabel('Water volume (m$^3$)')

ax += 1
axes[ax].plot(yr, lakeArea)
axes[ax].set_ylabel('Lake area (m$^2$)')

ax += 1
axes[ax].plot(yr, lakeVol)
axes[ax].set_ylabel('Lake volume (m$^3$)')

for ax in axes:
    ax.grid(True)
    ax.set_xlabel("Year")


plt.show()

