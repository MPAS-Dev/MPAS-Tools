#!/usr/bin/env python
'''
Script to fine tune ice-shelf geometry:
1. Smooth ice thickness across ice shelves.  The first row of cells adjacent to GL is left unmodified.
2. Adjust bathymetry to prevent spurious grounding from thickness adjustment by modifying 
bathymetry to maintain original water column thickness
3. Lower bathymetry beneath ice shelves and in open ocean by a uniform amount.  This reduces spurious 
grounding line advance from transient model dH/dt and geometry inconsistencies.  5 meters appears to be a good amount.

See beginning of script for options to adjust.

Note: changes are made to file in place so perform on a copy of your original file!

Matt Hoffman, 9/11/2022
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt


# ---- options to adjust (could become command line options ----

nSmoothing = 3 # number of rounds of smoothing to apply.  3 was the best match at 8km to 30 years of thickness evolution on Ross and FRIS.  Set to 0 to disable

elevationDrop = 5.0 # height in meters to drop non-groudned marine bed elevation

# ------------------------


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-f", dest="fileName", help="filename to modify", metavar="FILENAME")
options, args = parser.parse_args()

rhoi=910.0
rhosw=1028.0

f = Dataset(options.fileName, 'r+')

nCells = len(f.dimensions['nCells'])
mu = f.variables['muFriction'][0,:]
thickness = f.variables['thickness'][0,:]
bedTopography = f.variables['bedTopography'][0,:]
cOnC= f.variables['cellsOnCell'][:]
nEOnC = f.variables['nEdgesOnCell'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
lowerSurface = -rhoi/rhosw*thickness # only works for floating cells
WCT = lowerSurface - bedTopography

floatMask = ((thickness*910/1028+bedTopography)<0.0)*(thickness>0.0)
groundMask = ((thickness*910/1028+bedTopography)>0.0)*(thickness>0.0)

oGLmask = np.zeros((nCells,), 'i')
for c in range(nCells):
    if floatMask[c] == 1:
        for n in range(nEOnC[c]):
            if groundMask[cOnC[c,n]-1] == 1:
                oGLmask[c] = True
                break

oMgnMask = np.zeros((nCells,),'i')
for c in range(nCells):
    if floatMask[c] == 1:
        for n in range(nEOnC[c]):
            if thickness[cOnC[c,n]-1] == 0.0:
                oMgnMask[c] = True
                break

oMgnMaskOld = oMgnMask.copy()
for c in range(nCells):
    if floatMask[c] == 1:
        for n in range(nEOnC[c]):
            if oMgnMaskOld[cOnC[c,n]-1] == 1:
                oMgnMask[c] = True



mask = np.logical_and(np.logical_not(np.logical_or(oMgnMask, oGLmask)), floatMask) # mask of the floating area excluding the boundaries

#fig, axs = plt.subplots(1,1, figsize=(9,9), num=1)
#plt.scatter(xCell, yCell, s=2, c=mask*1+floatMask*1)
#plt.colorbar()
#plt.show()

# smooth over the mask

def smooth(H):
    Hold = H.copy()
    for c in range(nCells):
        if mask[c] == 1:
            Hsum = 0.0
            cnt = nEOnC[c]
            for n in range(cnt):
                Hsum += Hold[cOnC[c,n]-1]
            Hsum += Hold[c]
            cnt += 1
            H[c] = Hsum / float(cnt)
    return H

for i in range(nSmoothing):
  print(f"start smoothing {i}")
  thickness = smooth(thickness)
print("done")

f.variables['thickness'][0,:] = thickness

# Keep WCT
newLowerSurface = -rhoi/rhosw*thickness
for c in range(nCells):
    if floatMask[c] == 1:
        bedTopography[c] = newLowerSurface[c] - WCT[c]
        bedTopography[c] -= elevationDrop
f.variables['bedTopography'][0,:] = bedTopography

f.close()

print("saved")

fig, axs = plt.subplots(1,1, figsize=(9,9), num=1)
plt.scatter(xCell, yCell, s=2, c=thickness, vmin=200, vmax=1000.0, cmap='RdBu')
plt.colorbar()
plt.show()
