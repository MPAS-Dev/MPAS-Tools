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
dcEdge = f.variables['dcEdge'][:]
lowerSurface = -rhoi/rhosw*thickness # only works for floating cells
WCT = lowerSurface - bedTopography

floatMask = ((thickness*910/1028+bedTopography)<0.0)*(thickness>0.0)
#groundMask = ((thickness*910/1028+bedTopography)>0.0)*(thickness>0.0)


dCell = dcEdge.min()*2.5
dCell = 8000.0
print(f"using dCell={dCell}")

def smoothTrough(WCT, p1, p2, minWCT, maxWCT):
    WCTnew = WCT.copy()
    mask = (WCT<maxWCT) * (floatMask==1) * (xCell>(p1[0]-2*dCell)) * (xCell<(p2[0]+2*dCell)) * (yCell>(p1[1]-2*dCell)) * (yCell<(p2[1]+2*dCell))
    ind = np.nonzero(mask==1)[0]
    length = ( (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 )**0.5 # length along flowline
    m1 = (p2[1]-p1[1]) / (p2[0]-p1[0]) # xy slope of flowline 
    b1 = p1[1] - m1*p1[0] # y int of flowline
    print(f'# cells={len(ind)}, length={length}, m={m1}')
    m2 = -1.0/m1 # slope of all transverse lines
    for c in ind:
        # find width and center at this position
        b2 = yCell[c] - m2 * xCell[c]  # y int of transverse line
        dist2orthogline = np.absolute(m2*xCell + -1*yCell + b2) / (m2**2 + (-1)**2)**0.5
        ind2 = np.nonzero((mask==1) * (dist2orthogline<(1.5*dCell)))[0] # get the points within a certain distance of the orthog line
        dist2centerline = np.absolute(m1*xCell[ind2] + -1*yCell[ind2] + b1) / (m1**2 + (-1)**2)**0.5
        #dist2centerline = np.absolute( (p2[0]-p1[0]) * (p1[1] - yCell[ind2]) - (p1[0]-xCell[ind2])*(p2[1]-p1[1])) / ( (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)**0.5
        dist2centerline *= np.sign((m1*xCell[ind2]+b1) - yCell[ind2]) # make signed distance, assuming NE direction
        print(dist2centerline)
        width = dist2centerline.max()-dist2centerline.min()
        #width = dist2centerline.max()*2.0


        # find frac x along center line
        # first need intersection of longit. and ortho lines
        xint = (b2-b1)/(m1-(-1.0/m1))
        mag = (xint-p1[0]) / (p2[0]-p1[0]) * (maxWCT-minWCT) + minWCT


        widthOrgnIdx = ind2[np.argmin(dist2centerline)] #idx of westmost
        widthOrgn = (xCell[widthOrgnIdx], yCell[widthOrgnIdx]) # position of westmost
        widthEndIdx = ind2[np.argmax(dist2centerline)] #idx of westmost
        widthEnd = (xCell[widthEndIdx], yCell[widthEndIdx]) # position of westmost
        widthMidpt = ((widthOrgn[0]+widthEnd[0])/2.0, (widthOrgn[1]+widthEnd[1])/2.0)
        print(xCell[c], yCell[c], widthOrgn)
        #dist2Orgn = ((widthOrgn[0]-xCell[c])**2 + (widthOrgn[1]-yCell[c])**2)**0.5
        #fracWidth = (dist2Orgn  / (width) - 0.5) * 2.0 # [-1,1] range
        dist2Mid = ((widthMidpt[0]-xCell[c])**2 + (widthMidpt[1]-yCell[c])**2)**0.5
        fracWidth = dist2Mid  / (0.5*(width+2*dCell))  # [-1,1] range
        fracWidth = min(fracWidth, 0.95)
        print(f'found {len(ind2)} points along this orthogonal; width={width}; fracWidth={fracWidth}, mag={mag}')
        z = mag * (-1.0 * fracWidth**2 + 1.0)
        WCTnew[c] = max(z, WCT[c])  # don't make WCT thinner than the original
        print(f'old z={WCT[c]}, new z={z}')

    return WCTnew


maxWCT = 300.0
minWCT = 50.0
WCTnew = WCT.copy()
p1=[-1261073.6, 137133.0]; p2=[-1173862.4, 199332.15] # Rutford
WCTnew = smoothTrough(WCTnew, p1, p2, minWCT, maxWCT)
p1=[-1315459.94, 202346.75]; p2=[-1257332.3, 291441.3] # Carlson
WCTnew = smoothTrough(WCTnew, p1, p2, minWCT, maxWCT)
p1=[-1513787.6, 311996.4]; p2=[-1386639.4, 337107.2] # Evans
WCTnew = smoothTrough(WCTnew, p1, p2, minWCT, maxWCT)
p1=[-34186, 1983479]; p2=[-9070, 2047275] # Fimbul
WCTnew = smoothTrough(WCTnew, p1, p2, minWCT, maxWCT)
p1=[-1064519, 206317]; p2=[-1055596, 265431] # lil guy
WCTnew = smoothTrough(WCTnew, p1, p2, 40.0, 150.0)

f.variables['bedTopography'][0,:] = lowerSurface - WCTnew

f.close()

print("saved")
s=12
ind=(floatMask==1)
fig, axs = plt.subplots(1,2, figsize=(14,9), num=1, sharex=True, sharey=True)
plt.sca(axs[0])
plt.scatter(xCell[ind], yCell[ind], s=s, c=WCT[ind], vmin=0, vmax=100.0, cmap='RdBu')
plt.colorbar()
plt.sca(axs[1])
plt.scatter(xCell[ind], yCell[ind], s=s, c=WCTnew[ind], vmin=0, vmax=100.0, cmap='RdBu')
plt.plot(p1[0], p1[1], 'k*')
plt.plot(p2[0], p2[1], 'ko')
axs[0].set_aspect('equal', 'box')
axs[1].set_aspect('equal', 'box')
plt.colorbar()
plt.show()
