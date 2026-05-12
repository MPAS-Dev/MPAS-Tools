#!/usr/bin/env python
'''
Script to plot common time-series from one or more landice globalStats files.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import glob
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import textwrap

rhoi = 910.0
rhosw = 1028.
rhofw = 1000.

print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = ArgumentParser(description=__doc__)
parser.add_argument('files', nargs='*', default=['globalStats.nc'],
                    help='input filename(s)', metavar='FILENAME')
parser.add_argument("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="Gt", metavar="UNITS")
parser.add_argument("-c", dest="plotChange", help="plot time series as change from initial.  (not applied to GL flux or calving flux)  Without this option, the full magnitude of time series is used", action='store_true', default=False)
parser.add_argument("-s", dest="saveFile", help="file name to save png", default=None)
options = parser.parse_args()

# Expand any quoted glob patterns (e.g. "output*.nc")
expanded = []
for f in options.files:
    matches = glob.glob(f)
    if matches:
        expanded.extend(sorted(matches))
    else:
        expanded.append(f)

options.files = expanded
print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))


if options.units == "m3":
   massUnit = "m$^3$"
   fluxUnit = "m$^3$ yr$^{-1}$"
   scaleVol = 1.
   scaleFlux = 1.0 / rhoi
elif options.units == "kg":
   massUnit = "kg"
   fluxUnit = "kg yr$^{-1}$"
   scaleVol = rhoi
   scaleFlux = 1.0
elif options.units == "Gt":
   massUnit = "Gt"
   fluxUnit = "Gt yr$^{-1}$"
   scaleVol = rhoi / 1.0e12
   scaleFlux = 1.0 / 1.0e12
else:
   sys.exit("Unknown mass/volume units")
print("Using volume/mass units of: ", massUnit)

if options.plotChange:
    plotChangeStr = ' change'
else:
    plotChangeStr = ''

# ---
# Volume & area plots
figVol = plt.figure(1, figsize=(9, 12), facecolor='w')

nrow = 4
ncol = 2
ind = 1

axVol = figVol.add_subplot(nrow, ncol, ind)
plt.xlabel('Year')
plt.ylabel(f'volume{plotChangeStr} ({massUnit})')
plt.grid()
ind += 1
axX = axVol

axVAF = figVol.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'VAF{plotChangeStr} ({massUnit})')
plt.grid()
ind += 1

axVolGround = figVol.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'grounded volume{plotChangeStr} ({massUnit})')
plt.grid()
ind += 1

axVolFloat = figVol.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'floating volume{plotChangeStr} ({massUnit})')
plt.grid()
ind += 1

# area plots

axGrdArea = figVol.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'grounded area{plotChangeStr} (km$^2$)')
plt.grid()
ind += 1

axFltArea = figVol.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'floating area{plotChangeStr} (km$^2$)')
plt.grid()
ind += 1

axTotArea = figVol.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'total area{plotChangeStr} (km$^2$)')
plt.grid()
ind += 1

axLeg = figVol.add_subplot(nrow, ncol, ind)
axLeg.axis('off')

# ---
# Flux plots
figFlux = plt.figure(2, figsize=(9, 12), facecolor='w')

nrow = 5
ncol = 2
ind = 1

axTotSMBFlux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'total SMB flux ({fluxUnit})')
plt.grid()
ind += 1

axTotBMBFlux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'total BMB flux ({fluxUnit})')
plt.grid()
ind += 1

axGrdSMBFlux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'grounded SMB flux ({fluxUnit})')
plt.grid()
ind += 1

axGrdBMBFlux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'grounded BMB flux ({fluxUnit})')
plt.grid()
ind += 1

axFltSMBFlux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'floating SMB flux ({fluxUnit})')
plt.grid()
ind += 1

axFltBMBFlux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'floating BMB flux ({fluxUnit})')
plt.grid()
ind += 1

axGLflux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'GL flux ({fluxUnit})')
plt.grid()
ind += 1

axGLMigflux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'GL migration flux ({fluxUnit})')
plt.grid()
ind += 1

axCalvFlux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'calving flux ({fluxUnit})')
plt.grid()
ind += 1

axFacemeltFlux = figFlux.add_subplot(nrow, ncol, ind, sharex=axX)
plt.xlabel('Year')
plt.ylabel(f'facemelt flux ({fluxUnit})')
plt.grid()
ind += 1


def VAF2seaLevel(vol):
    """
    This function accounts for the fact that ice, when melted,
    takes on freshwater density.  This density correction is only
    applied to the volume above flotation, when in reality, it should
    be applied for all melted ice.
    """
    return vol / scaleVol / 3.62e14 * rhoi / rhofw * 1000.

def seaLevel2VAF(vol):
    return vol * scaleVol * 3.62e14 * rhofw / rhoi / 1000.

def addSeaLevAx(axName):
    seaLevAx = axName.secondary_yaxis('right', functions=(VAF2seaLevel, seaLevel2VAF))
    seaLevAx.set_ylabel('Sea-level\nequivalent (mm)')

def custom_legend():
    handles, labels = axVol.get_legend_handles_labels()
    leg = axLeg.legend(handles, labels, loc='center', prop={'size': 7})
    # add line breaks
    max_width = 65
    leg_texts = leg.get_texts()
    for ind in range(len(leg_texts)):
        text = leg_texts[ind].get_text()
        if len(text) > max_width:
            wrapped_text = '\n'.join(textwrap.wrap(text, width=max_width, break_on_hyphens=False, subsequent_indent='  '))
            leg_texts[ind].set(text=wrapped_text)

def plotStat(fname):
    print("Reading and plotting file: {}".format(fname))

    name = fname

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    yr = yr-yr[0]
    dt = f.variables['deltat'][:]/3.15e7
    print(f"Year range: {yr.min()}, {yr.max()}")

    # volume stats

    vol = f.variables['totalIceVolume'][:] * scaleVol
    if options.plotChange:
        vol = vol - vol[0]
    axVol.plot(yr, vol, label=name)

    VAF = f.variables['volumeAboveFloatation'][:] * scaleVol
    if options.plotChange:
        VAF = VAF - VAF[0]
    axVAF.plot(yr, VAF, label=name)
    addSeaLevAx(axVAF)

    volGround = f.variables['groundedIceVolume'][:] * scaleVol
    if options.plotChange:
        volGround = volGround - volGround[0]
    axVolGround.plot(yr, volGround, label=name)

    volFloat = f.variables['floatingIceVolume'][:] * scaleVol
    if options.plotChange:
        volFloat = volFloat - volFloat[0]
    axVolFloat.plot(yr, volFloat, label=name)

    # area stats

    areaGrd = f.variables['groundedIceArea'][:] / 1000.0**2
    if options.plotChange:
        areaGrd = areaGrd - areaGrd[0]
    axGrdArea.plot(yr, areaGrd, label=name)

    areaFlt = f.variables['floatingIceArea'][:] / 1000.0**2
    if options.plotChange:
        areaFlt = areaFlt - areaFlt[0]
    axFltArea.plot(yr, areaFlt, label=name)

    areaTot = f.variables['totalIceArea'][:] / 1000.0**2
    if options.plotChange:
        areaTot = areaTot - areaTot[0]
    axTotArea.plot(yr, areaTot, label=name)

    # flux stats

    totSMB = f.variables['totalSfcMassBal'][:] * scaleFlux
    axTotSMBFlux.plot(yr, totSMB, label=name)

    totBMB = f.variables['totalBasalMassBal'][:] * scaleFlux
    axTotBMBFlux.plot(yr, totBMB, label=name)

    grdSMB = f.variables['totalGroundedSfcMassBal'][:] * scaleFlux
    axGrdSMBFlux.plot(yr, grdSMB, label=name)

    grdBMB = f.variables['totalGroundedBasalMassBal'][:] * scaleFlux
    axGrdBMBFlux.plot(yr, grdBMB, label=name)

    fltSMB = f.variables['totalFloatingSfcMassBal'][:] * scaleFlux
    axFltSMBFlux.plot(yr, fltSMB, label=name)

    fltBMB = f.variables['totalFloatingBasalMassBal'][:] * scaleFlux
    axFltBMBFlux.plot(yr, fltBMB, label=name)

    GLflux = f.variables['groundingLineFlux'][:] * scaleFlux
    axGLflux.plot(yr, GLflux, label=name)

    GLMigflux = f.variables['groundingLineMigrationFlux'][:] * scaleFlux
    axGLMigflux.plot(yr, GLMigflux, label=name)

    calvFlux = f.variables['totalCalvingFlux'][:] * scaleFlux
    axCalvFlux.plot(yr, calvFlux, label=name)

    fmFlux = f.variables['totalFaceMeltingFlux'][:] * scaleFlux
    axFacemeltFlux.plot(yr, fmFlux, label=name)

    f.close()


for fname in options.files:
    plotStat(fname)

custom_legend()

print("Generating plot.")
figVol.tight_layout()
figFlux.tight_layout()
if options.saveFile is not None:
    figVol.savefig(f"{options.saveFile}_vol.png", dpi=400, bbox_inches='tight')
    figFlux.savefig(f"{options.saveFile}_flux.png", dpi=400, bbox_inches='tight')
plt.show()

