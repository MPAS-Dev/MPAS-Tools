#!/usr/bin/env python
'''
Script to plot mass balance time-series from landice globalStats files.
Currently only assesses grounded ice sheet mass balance.
'''

import sys
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt

rhoi = 910.0


print "** Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser(description=__doc__)
parser.add_option("-f", dest="fileName", help="input filename", default="globalStats.nc", metavar="FILENAME")
options, args = parser.parse_args()

print "Using ice density of {} kg/m3 if required for unit conversions".format(rhoi)

print "Mass balance will be inaccurate if not writing stats on every timestep."

print "Reading and plotting file: " + options.fileName
f = Dataset(options.fileName,'r')
yr = f.variables['daysSinceStart'][:]/365.0
dyr = np.zeros(yr.shape)
dyr[1:] = yr[1:]-yr[:-1]

# ---- Figure: mass change rate ---
fig = plt.figure(1, figsize=(6, 5), facecolor='w')

ax = fig.add_subplot(1, 1, 1)
plt.xlabel('Year')
plt.ylabel('Grounded mass change (Gt/yr)')
plt.grid()


# get fields and plot stuff
volGround = f.variables['groundedIceVolume'][:]
volGround = volGround * rhoi / 1.0e12
dvolGround = np.zeros(yr.shape)
dvolGround[1:] = (volGround[1:]-volGround[:-1]) / dyr[1:]
ax.plot(yr, dvolGround, 'k', linewidth=3, label="grounded mass change")

SMBg = f.variables['totalGroundedSfcMassBal'][:] / 1.0e12
ax.plot(yr, SMBg, label="grounded SMB")

BMBg = f.variables['totalGroundedBasalMassBal'][:] / 1.0e12
ax.plot(yr, BMBg, label="grounded BMB")

GLflux = f.variables['groundingLineFlux'][:]/1.0e12
ax.plot(yr, -GLflux, label="GL flux")

GLMigrationflux = f.variables['groundingLineMigrationFlux'][:]/1.0e12
ax.plot(yr, -GLMigrationflux, label="GL migration flux")

ax.plot(yr, SMBg+BMBg, "--", label="SMB+BMB")
ax.plot(yr, SMBg+BMBg-GLflux, "--", label="SMB+BMB+GLf")
ax.plot(yr, SMBg+BMBg-GLflux-GLMigrationflux, "--", label="SMB+BMB+GLf+GLMf")

plt.legend(loc='best', prop={'size': 6})

# --- Figure: cumulative mass change ---
fig2 = plt.figure(2, figsize=(6,5), facecolor='w')
ax = fig.add_subplot(1, 1, 1)
plt.xlabel('Year')
plt.ylabel('Cumulative grounded mass change (Gt)')
plt.grid()

plt.plot(yr, volGround - volGround[0], 'k', linewidth=3, label='grounded mass')
plt.plot(yr, SMBg.cumsum(), label="grounded SMB")
plt.plot(yr, BMBg.cumsum(), label="grounded BMB")
plt.plot(yr, GLflux.cumsum(), label="GL flux")
plt.plot(yr, GLMigrationflux.cumsum(), label="GL Migration flux")
plt.plot(yr, ( (SMBg+BMBg) * dyr).cumsum(), "--", label='SMB+BMB')
plt.plot(yr, ( (SMBg+BMBg-GLflux) * dyr).cumsum(), "--", label='SMB+BMB+GLflux')
plt.plot(yr, ( (SMBg+BMBg-GLflux-GLMigrationflux) * dyr).cumsum(), "--", label='SMB+BMB+GLf+GLMf')
#plt.plot(yr, dareag.cumsum(), 'k', label="mass change due to change in grounded area")

plt.legend(loc='best', prop={'size': 6})

print "Plotting complete."

plt.show()
f.close()
