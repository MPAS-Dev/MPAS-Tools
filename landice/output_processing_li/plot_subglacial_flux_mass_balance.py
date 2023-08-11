#!/usr/bin/env python
'''
Script to plot subglacial hydrology time-series from a landice globalStats files.
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

totalSubglacialWaterVolume = dataset.variables['totalSubglacialWaterVolume'][:]
melt = dataset.variables['totalBasalMeltInput'][:]
fluxMarine = dataset.variables['totalGLMeltFlux'][:]
fluxLand = dataset.variables['totalTerrestrialMeltFlux'][:]


chnlFluxMarine = dataset.variables['totalChannelGLMeltFlux'][:]
chnlFluxLand = dataset.variables['totalChannelTerrestrialMeltFlux'][:]
chnlMelt = dataset.variables['totalChannelMelt'][:]

deltat = dataset.variables['deltat'][:]
yr = dataset.variables['daysSinceStart'][:] / 365.0

subglacialWaterMassRate = np.zeros((len(melt),))
rhow = 1000.0

for i in range(len(totalSubglacialWaterVolume) - 1):
    subglacialWaterMassRate[i] = ((totalSubglacialWaterVolume[i+1] - totalSubglacialWaterVolume[i]) * rhow / deltat[i])



#print(subglacialWaterMassRate)

# input
plt.plot(yr, melt, 'r:', label='basal melt')
plt.plot(yr, chnlMelt, 'r--', label='channel melt')
total_melt = melt + chnlMelt
plt.plot(yr, total_melt, 'r-', label='total melt')

# output
plt.plot(yr, fluxMarine+fluxLand, 'b:', label='sheet outflux')
plt.plot(yr, chnlFluxMarine+chnlFluxLand, 'b--', label='chnl outflux')
total_outflux = fluxMarine+fluxLand + chnlFluxMarine+chnlFluxLand
plt.plot(yr, total_outflux, 'b-', label='total outflux')

plt.plot(yr, subglacialWaterMassRate, 'g-', label='dV/dt')

plt.plot(yr, total_melt - total_outflux, 'k', label='I-O')

plt.legend()
plt.xlabel('Year')
plt.ylabel('Mass flux (kg/s)')
plt.show()

