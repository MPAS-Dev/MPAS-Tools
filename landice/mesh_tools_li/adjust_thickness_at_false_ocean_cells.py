#!/usr/bin/env python
'''
This script changes the designation of isolated cells inland from the grounding line that are wrongfully defined as 
floating ice (and thus part of the ocean) when using the ocean density to define grounded/floating cells, as done in li_mask_is_grounded_ice.
Ocean and floating ice cells are flood filled from the edges of the domain, and any floating ice cells not in contact with the flood fill are
identified. Thickness of these cells is then manually altered to enforce the designation of grounded ice.
'''
import mpas_tools
import numpy as np
import xarray as xr
from compass.landice.mesh import mpas_flood_fill
from optparse import OptionParser
import subprocess

print("** Gathering information ...")
parser = OptionParser()
parser.add_option("-f", "--file", dest="file", metavar="FILE")
options, args = parser.parse_args()

f = xr.open_dataset(options.file, decode_times=False, decode_cf=False)
cellsOnCell= f['cellsOnCell'][:,:].data
nEdgesOnCell = f['nEdgesOnCell'][:].data
thickness = f['thickness'][0,:].data
bedTopography = f['bedTopography'][0,:].data
cellsOnEdge = f['cellsOnEdge'][:,:].data

groundedIceMask = ((thickness*910/1028+bedTopography)>0.0)*(thickness>0.0)
floatingIceMask = ((thickness*910/1028+bedTopography)<=0.0)*(thickness>0.0)
oceanMask = (thickness==0.0)*(bedTopography<0.0)
landMask = (thickness==0.0)*(bedTopography>=0.0)

seedMask = np.zeros((len(nEdgesOnCell),), 'float64')
growMask = floatingIceMask + oceanMask

ind = np.where(oceanMask==1)[0]
for i in ind:
    for ii in range(nEdgesOnCell[i]):
        if (cellsOnCell[i,ii] == 0):
            seedMask[i] = 1

print("**Flood Filling ...")

keepMask = mpas_flood_fill(seedMask, growMask, cellsOnCell, nEdgesOnCell)

ind = np.where(floatingIceMask==1)[0]
for i in ind:
    if (keepMask[i] == 0):
        thickness[i] = -bedTopography[i] * 1028/910 + 1e-10 #thickness necessary to achieve grounded ice (small margin past exact flotation)

print("**Saving ...")

seedMask = seedMask.reshape(1,len(seedMask))
growMask = growMask.reshape(1,len(growMask))
keepMask = keepMask.reshape(1,len(keepMask))
thickness = thickness.reshape(1,len(thickness))

gm = xr.DataArray(growMask.astype('float64'),dims=('Time','nCells'))
f['growMask'] = gm

sm = xr.DataArray(seedMask.astype('float64'),dims=('Time','nCells'))
f['seedMask'] = sm

km = xr.DataArray(keepMask.astype('float64'),dims=('Time','nCells'))
f['keepMask'] = km

thk = xr.DataArray(thickness.astype('float64'),dims=('Time','nCells'))
f['thickness'] = thk

f.to_netcdf("modifiedThicknessDomain.nc")
f.close()

subprocess.run(["ncatted", "-a", "_FillValue,,d,,", "modifiedThicknessDomain.nc"])
