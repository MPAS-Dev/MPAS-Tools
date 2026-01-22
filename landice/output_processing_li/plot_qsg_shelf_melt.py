#!/usr/bin/env python

import numpy as np
import xarray as xr
from  matplotlib import pyplot as plt
from optparse import OptionParser
from pyproj import Transformer, transform, CRS
from scipy.interpolate import griddata
import cmocean

parser = OptionParser()
parser.add_option("-m", "--oceanMesh", dest="oceanMesh", help="Filename with Mpas-Ocean mesh", metavar="FILENAME")
parser.add_option("-o", "--oceanOutput", dest="oceanOutput", help="Filename with Mpas-Ocean output", metavar="FILENAME")
parser.add_option("-r", "--region", dest="region", help="Name of region to plot. Options are 'Ross', 'Peninsula', 'Amundsen', 'Filchner-Ronne', 'Amery', 'AIS'", default="AIS")
options, args = parser.parse_args()

#<TO DO>: Add Qsg throughout

# load variables
OM = xr.open_dataset(options.oceanMesh,decode_timedelta=False)
DS = xr.open_dataset(options.oceanOutput,decode_timedelta=False)

lonCell = OM['lonCell'].data
latCell = OM['latCell'].data
melt = DS['timeMonthly_avg_landIceFreshwaterFlux'].data

#isolate antarctic
latCellDeg = latCell * 180 / np.pi
lonCellDeg = lonCell * 180 / np.pi

ind = np.where(latCellDeg <= -65)
latCellAnt = latCellDeg[ind]
lonCellAnt = lonCellDeg[ind]
meltAnt = melt[0,ind]

#convert to polar stereographic
projections = dict()
projections['ais-bedmap2'] = '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84' # Note: BEDMAP2 elevations use EIGEN-GL04C geoid
projections['latlon'] = '+proj=longlat +ellps=WGS84'

crs_out = CRS.from_proj4(projections['ais-bedmap2'])
crs_in = CRS.from_proj4(projections['latlon'])

print("crs_out: {}".format(crs_out))
print("crs_in: {}".format(crs_in))

# define a transformer
t = Transformer.from_crs(crs_in,crs_out)

xCellAnt,yCellAnt = t.transform(lonCellAnt[:], latCellAnt[:], radians = False)

# Set up geographic windows
print("Setting up " + options.region + " Region")
if (options.region == 'Ross'):
    xlimit = (-10e5, 8e5)
    ylimit = (-2.2e6, -0.4e6)
elif (options.region == 'Peninsula'):
    xlimit = (-2.7e6, -1.3e6)
    ylimit = (0, 18e5)
elif (options.region == 'Amundsen'):
    xlimit = (-2.1e6, -0.4e6)
    ylimit = (-16e5, 0)
elif (options.region == 'Filchner-Ronne'):
    xlimit = (-17e5, -4e5)
    ylimit = (1.3e5, 15e5)
elif (options.region == 'Amery'):
    xlimit = (1.4e6, 2.8e6)
    ylimit = (-2e5, 14e5)
else : #AIS
    xlimit = (-2.75e6, 2.75e6)
    ylimit = (-2.25e6, 2.25e6)

#isolate x/y indices within region
ind = np.array(np.where((xCellAnt >= np.min(xlimit)) & (xCellAnt <= np.max(xlimit)) & (yCellAnt >= np.min(ylimit)) & (yCellAnt <= np.max(ylimit))))
#set up grid to interpolate onto - better for visualization purposes
res = 15*1000 #15 km grid resolution
Xvec = np.arange(np.min(xCellAnt[ind]), np.max(xCellAnt[ind]) + res, res)
Yvec = np.arange(np.min(yCellAnt[ind]), np.max(yCellAnt[ind]) + res, res)

Xgrd, Ygrd = np.meshgrid(Xvec,Yvec)
Xgrd = np.transpose(Xgrd)
Ygrd = np.transpose(Ygrd)

xp = np.array(xCellAnt[ind])
yp = np.array(yCellAnt[ind])

values = melt[0,ind]

Mgrd = np.full_like(Xgrd, np.nan)

#bin average onto grid
for i in range(len(Xvec)-1):
    for ii in range(len(Yvec)-1):
        ind = np.array(np.where((xp >= Xvec[i]) & (xp <= Xvec[i+1]) & (yp >= Yvec[ii]) & (yp <= Yvec[ii+1])))
        if ind.size != 0:
            Mgrd[i,ii] = np.mean(values[0,ind])
        else:
            Mgrd[i,ii] = np.nan

Mgrd = np.squeeze(Mgrd)
Mgrd = Mgrd*3.1536e7/1000 #convert to m/yr
# Create plots
plt.pcolormesh(Xgrd,Ygrd,Mgrd,shading='gouraud',vmin=0,vmax=0.5)
plt.colorbar()
plt.show()
plt.savefig('test.png')
