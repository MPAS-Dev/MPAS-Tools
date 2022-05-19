#!/usr/bin/env python
'''
Script to plot mass balance time-series from landice globalStats files.
Currently only assesses grounded ice sheet mass balance.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt

rhoi = 910.0


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-f", dest="fileName", help="input filename", default="globalStats.nc", metavar="FILENAME")
options, args = parser.parse_args()

print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))

print("Mass balance will be inaccurate if not writing stats on every timestep.")

print("Reading and plotting file: {}".format(options.fileName))
f = Dataset(options.fileName,'r')
yr = f.variables['daysSinceStart'][:]/365.0
dyr = np.zeros(yr.shape)
dyr[1:] = yr[1:]-yr[:-1]

# ---- Figure: mass change rate ---
fig, ax = plt.subplots(3,3, figsize=(15,12))
#=====================
# 1. Total mass budget
#=====================
# get fields and plot stuff
vol = f.variables['totalIceVolume'][:]
vol = vol * rhoi / 1.0e12
dvol = np.zeros(yr.shape)
dvol[1:] = (vol[1:]-vol[:-1]) / dyr[1:]

# --- time-step-wise mass change ---
ax[0,0].plot(yr, dvol, 'k', linewidth=3, label="mass change")

SMB = f.variables['totalSfcMassBal'][:] / 1.0e12
ax[0,0].plot(yr, SMB, label="SMB")

BMB = f.variables['totalBasalMassBal'][:] / 1.0e12
ax[0,0].plot(yr, BMB, label="BMB")

calv = -f.variables['totalCalvingFlux'][:] / 1.0e12
ax[0,0].plot(yr, calv, label='calv')

FMF = -f.variables['totalFaceMeltingFlux'][:] / 1.0e12
ax[0,0].plot(yr, FMF, label='facemelt')

# Grounding line flux and grounding line migration flux
# do not come into the global mass budget
tot = SMB+BMB+calv+FMF
ax[0,0].plot(yr, tot, "--", label="total")

ax[0,0].legend(loc='best', prop={'size': 6})
ax[0,0].set_ylabel('Mass change (Gt)')
# --- cumulative mass change ---

ax[1,0].plot(yr, vol - vol[0], 'k', linewidth=3, label='total mass change')
ax[1,0].plot(yr, (SMB*dyr).cumsum(), label="SMB")
ax[1,0].plot(yr, (BMB*dyr).cumsum(), label="BMB")

ax[1,0].plot(yr, (calv*dyr).cumsum(), label='calving')
ax[1,0].plot(yr, (FMF*dyr).cumsum(), label='facemelt')
ax[1,0].plot(yr, (tot * dyr).cumsum(), "--", label='total budget')

ax[1,0].legend(loc='best', prop={'size': 6})
ax[1,0].set_ylabel('Cumulative mass change (Gt)')

ax[2,0].semilogy(yr, np.abs( (tot - dvol ) / dvol ) )
ax[2,0].set_ylabel('fractional error', color='tab:blue')
absErrAxTot = ax[2,0].twinx()
absErrAxTot.plot(yr, tot - dvol, color='tab:orange')
#========================
# 2. Grounded mass budget
#========================
# get fields and plot stuff
volGround = f.variables['groundedIceVolume'][:]
volGround = volGround * rhoi / 1.0e12
dvolGround = np.zeros(yr.shape)
dvolGround[1:] = (volGround[1:]-volGround[:-1]) / dyr[1:]
ax[0,1].plot(yr, dvolGround, 'k', linewidth=3, label="mass change")

# --- time-step-wise mass change ---
SMBg = f.variables['totalGroundedSfcMassBal'][:] / 1.0e12
ax[0,1].plot(yr, SMBg, label="SMB")

BMBg = f.variables['totalGroundedBasalMassBal'][:] / 1.0e12
ax[0,1].plot(yr, BMBg, label="BMB")

calvg = -f.variables['totalGroundedCalvingFlux'][:] / 1.0e12
ax[0,1].plot(yr, calvg, label='calv')

FMFg = -f.variables['totalGroundedFaceMeltingFlux'][:] / 1.0e12
ax[0,1].plot(yr, FMFg, label='facemelt')

GLflux = -f.variables['groundingLineFlux'][:]/1.0e12
ax[0,1].plot(yr, GLflux, label="GL flux")

GLMigrationflux = -f.variables['groundingLineMigrationFlux'][:]/1.0e12
ax[0,1].plot(yr, GLMigrationflux, label="GL migration flux")

grndTot = SMBg+calvg+BMBg+FMFg+GLflux+GLMigrationflux
ax[0,1].plot(yr, grndTot, "--", label="total")

ax[0,1].legend(loc='best', prop={'size': 6})

# --- Figure: cumulative mass change ---

ax[1,1].plot(yr, volGround - volGround[0], 'k', linewidth=3, label='total mass change')
ax[1,1].plot(yr, (SMBg*dyr).cumsum(), label="SMB")
ax[1,1].plot(yr, (BMBg*dyr).cumsum(), label="BMB")

ax[1,1].plot(yr, (calvg*dyr).cumsum(), label='calving')
ax[1,1].plot(yr, (FMFg*dyr).cumsum(), label='facemelt')
ax[1,1].plot(yr, (GLflux*dyr).cumsum(), label="GL flux")
ax[1,1].plot(yr, (GLMigrationflux*dyr).cumsum(), label="GL Migration flux")
ax[1,1].plot(yr, ( (grndTot) * dyr).cumsum(), "--", label='total budget')

ax[1,1].legend(loc='best', prop={'size': 6})
ax[2,1].semilogy(yr, np.abs( (grndTot - dvolGround)  / dvolGround ) )

absErrAxGrnd = ax[2,1].twinx()
absErrAxGrnd.plot(yr, grndTot - dvolGround, color='tab:orange')
#========================
# 3. Floating mass budget
#========================

# get fields and plot stuff
volFloat = f.variables['floatingIceVolume'][:]
volFloat = volFloat * rhoi / 1.0e12
dvolFloat = np.zeros(yr.shape)
dvolFloat[1:] = (volFloat[1:]-volFloat[:-1]) / dyr[1:]
ax[0,2].plot(yr, dvolFloat, 'k', linewidth=3, label="mass change")

SMBf = f.variables['totalFloatingSfcMassBal'][:] / 1.0e12
ax[0,2].plot(yr, SMBf, label="SMB")

BMBf = f.variables['totalFloatingBasalMassBal'][:] / 1.0e12
ax[0,2].plot(yr, BMBf, label="BMB")

calvf = -f.variables['totalFloatingCalvingFlux'][:] / 1.0e12
ax[0,2].plot(yr, calvf, label='calv')

FMFf = -f.variables['totalFloatingFaceMeltingFlux'][:] / 1.0e12
ax[0,2].plot(yr, FMFf, label='facemelt')

ax[0,2].plot(yr, -GLflux, label="GL flux")

ax[0,2].plot(yr, -GLMigrationflux, label="GL migration flux")

fltTot = SMBf+BMBf+FMFf+calvf-GLflux-GLMigrationflux
ax[0,2].plot(yr, fltTot, "--", label="total")

ax[0,2].legend(loc='best', prop={'size': 6})

# --- cumulative mass change ---
ax[1,2].plot(yr, volFloat - volFloat[0], 'k', linewidth=3, label='total mass change')
ax[1,2].plot(yr, (SMBf*dyr).cumsum(), label="SMB")
ax[1,2].plot(yr, (BMBf*dyr).cumsum(), label="BMB")

ax[1,2].plot(yr, (calvf*dyr).cumsum(), label='calving')
ax[1,2].plot(yr, (FMFf*dyr).cumsum(), label='facemelt')
ax[1,2].plot(yr, (-GLflux*dyr).cumsum(), label="GL flux")
ax[1,2].plot(yr, (-GLMigrationflux*dyr).cumsum(), label="GL Migration flux")
ax[1,2].plot(yr, (fltTot * dyr).cumsum(), "--", label='total budget')

ax[1,2].legend(loc='best', prop={'size': 6})

ax[2,2].semilogy(yr, np.abs( (fltTot - dvolFloat) / dvolFloat ) )

absErrAxFlt = ax[2,2].twinx()
absErrAxFlt.plot(yr, fltTot - dvolFloat, color='tab:orange')
absErrAxFlt.set_ylabel('Absolute error (Gt; budget - true)', color='tab:orange')

for plotAx in ax.ravel():
    plotAx.grid('on')
ax[2,0].set_xlabel('yr')
ax[2,1].set_xlabel('yr')
ax[2,2].set_xlabel('yr')

ax[0,0].set_title('Total budget')
ax[0,1].set_title('Grounded budget')
ax[0,2].set_title('Floating budget')

fig.subplots_adjust(wspace=0.5)
plt.show()
f.close()
