#!/usr/bin/env python
'''
Uses output from MALI's thermal solver to determine where regions of the bed are frozen or thawed and masks off frozen cells from the hydro domain.
Only thawed cells in hydraulic contact with the grounding line are incorporated into the hydro domain. All other cells are set to zero ice thickness
and are surrounded by no-flow conditions with the waterFluxMask. An additional argument, "-u", will count any frozen ice as thawed if basal sliding speeds
exceed a threshold value (default is 25 m/yr).
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
parser.add_option("-u", "--UbThresh", dest="UbThresh", type="float", default=25, help="basal sliding threshold used to redefine frozen ice as thawed where sliding speed is above UbThresh (in m/yr)")
options, args = parser.parse_args()

f = xr.open_dataset(options.file, decode_times=False, decode_cf=False)
groundedBasalMassBal = f['groundedBasalMassBal'][0,:].data
cellsOnCell= f['cellsOnCell'][:,:].data
nEdgesOnCell = f['nEdgesOnCell'][:].data
thickness = f['thickness'][0,:].data
bedTopography = f['bedTopography'][0,:].data
edgesOnCell = f['edgesOnCell'][:,:].data
cellsOnEdge = f['cellsOnEdge'][:,:].data
yEdge = f['yEdge'][:].data
nVertInterfaces = f.sizes['nVertInterfaces']
uReconstructX = f['uReconstructX'][0,:,nVertInterfaces-1].data
uReconstructY = f['uReconstructY'][0,:,nVertInterfaces-1].data

groundedMask = ((thickness*910/1028+bedTopography)>0.0)*(thickness>0.0)
floatingMask = ((thickness*910/1028+bedTopography)<=0.0)*(thickness>0.0)
oceanMask = (thickness==0.0)*(bedTopography<0.0)
landMask = (thickness==0.0)*(bedTopography>=0.0)

groundedMask = groundedMask.reshape((len(groundedMask),1))
floatingMask = floatingMask.reshape(len(floatingMask),1)
oceanMask = oceanMask.reshape(len(oceanMask),1)
landMask = landMask.reshape(len(landMask),1)

seedMask = np.zeros((len(nEdgesOnCell),1), 'float64')

print("**Defining Seed and Grow Masks ...")
ind = np.where(groundedMask==1)[0]
for i in ind:
    #identify grounded cells just inland of grounding line
    for ii in range(nEdgesOnCell[i]):
       if ((floatingMask[cellsOnCell[i,ii]-1] == 1) | (oceanMask[cellsOnCell[i,ii]-1] == 1)):
          seedMask[i] = 1

#identify grounded, thawed ice
basalSlidingSpeed = np.sqrt(uReconstructX**2 + uReconstructY**2) * 3.15e7 #convert to m/yr

growMask = np.logical_or(groundedBasalMassBal < 0.0, basalSlidingSpeed > options.UbThresh)
growMask = growMask.reshape(len(growMask),1)

print("**Flood Filling ...")
keepMask = mpas_flood_fill(seedMask, growMask, cellsOnCell, nEdgesOnCell)

print("Defining waterFluxMask ...")

#zero out ice thickness where frozen ice or no thawed ice in contact with grounding line
thickness = thickness.reshape(len(thickness),1)
thickness[keepMask!=1] = 0

#create zero flux mask on edges of inactive domain
waterFluxMask = np.zeros((len(yEdge),1),'int32')
for i in range(len(yEdge)):
    cell1 = cellsOnEdge[i,0]-1
    cell2 = cellsOnEdge[i,1]-1
    if ((keepMask[cell1] == 1 and keepMask[cell2] == 0 and groundedMask[cell1] == 1 and groundedMask[cell2] == 1) \
        or (keepMask[cell2] == 1 and keepMask[cell1] == 0 and groundedMask[cell1] == 1 and groundedMask[cell2] == 1)):
        waterFluxMask[i] = 2

print("Saving ....")
try:
   f = f.drop_vars(['xtime'])
finally:
   print("No xtime variable to delete")

try:
   f = f.drop_vars(['simulationStartTime'])
finally:
   print("No simulationStartTime variable to delete")

try:
   f = f.drop_vars(['forcingTimeStamp'])
finally:
   ("No forcingTimeStamp variable to delete")

basalMeltInput = -1.0 * np.minimum(0.0, groundedBasalMassBal)
basalMeltInput = basalMeltInput.reshape(len(basalMeltInput), 1)
bmi = xr.DataArray(basalMeltInput.T.astype('float64'),dims=('Time','nCells'))
f['basalMeltInput'] = bmi

wfm = xr.DataArray(waterFluxMask.T.astype('int32'),dims=('Time','nEdges'))
f['waterFluxMask'] = wfm

thk = xr.DataArray(thickness.T.astype('float64'),dims=('Time','nCells'))
f['thickness'] = thk

f.to_netcdf("finalMaskedFile.nc")
f.close()

subprocess.run(["ncatted", "-a", "_FillValue,,d,,", "finalMaskedFile.nc"])
