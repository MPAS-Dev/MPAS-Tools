#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot bed and surface topography, surface speed, and
(optionally) temperature from MPAS netCDF along a transect.
@author: trevorhillebrand
"""

import numpy as np
import csv
from netCDF4 import Dataset
from optparse import OptionParser
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

parser = OptionParser(description='Plot transect from MPAS netCDF')
parser.add_option("-d", "--data", dest="data_file",
                  help="the MPAS netCDF file")
parser.add_option("-c", "--coords", dest="coords_file", default=None,
                  help="csv file defining transect. x coordinates \
                        in first column, y coordinates in second column. \
                        No header.")
parser.add_option('-t', dest="times", default="-1",
                  help="integer time levels at which to plot \
                        (int separated by commas; no spaces)")
parser.add_option('-x', dest='x_coords', default=None,
                  help='List of x coordinates of transect if not \
                        using a csv file. Comma-separated, no spaces.')
parser.add_option('-y', dest='y_coords', default=None,
                  help='List of y coordinates of transect if not \
                        using a csv file. Comma-separated, no spaces.')
parser.add_option('-s', dest='save_filename', default=None,
                  help='Path to save .png to, if desired.')

parser.add_option("--temperature", dest="interp_temp",
                  action="store_true", help="interpolate temperature")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

times_list = [i for i in options.times.split(',')]  # list of string times for plotting
times = [int(i) for i in options.times.split(',')]  # list of integer time indices

dataset = Dataset(options.data_file, 'r')
dataset.set_always_mask(False)
xCell = dataset.variables["xCell"][:]
yCell = dataset.variables["yCell"][:]
nCells = dataset.dimensions['nCells'].size
areaCell = dataset.variables["areaCell"][:]
layerThicknessFractions = dataset.variables["layerThicknessFractions"][:]
nVertLevels = dataset.dimensions['nVertLevels'].size
if "daysSinceStart" in dataset.variables.keys():
    use_yrs = True
    yrs = dataset.variables["daysSinceStart"][:] / 365.
    times_list = [f'{yrs[i]:.1f}' for i in times]
else:
    use_yrs = False

# Replace -1 time index with last forward index of time array.
# It is unclear why these need to be in separate if-statements, but they do.
if -1 in times:
    times[times.index(-1)] = int(dataset.dimensions['Time'].size - 1)
if '-1' in times_list:
    times_list[times_list.index('-1')] = str(dataset.dimensions['Time'].size - 1)

# Cannot plot temperature for more than one time index.
if options.interp_temp and (len(times) > 1):
    print('Cannot plot temperature for more than one time index.' +
          ' Skipping temperature interpolation and plotting.')
    options.interp_temp = False

li_mask_ValueDynamicIce = 2
cellMask = dataset.variables['cellMask'][times,:]
cellMask_dynamicIce = (cellMask & li_mask_ValueDynamicIce) // li_mask_ValueDynamicIce
# only take thickness of dynamic ice
thk = dataset.variables["thickness"][times,:] * cellMask_dynamicIce
plot_speed = True
# Include speed on non-dynamic ice to avoid interpolation artifacts.
if "surfaceSpeed" in dataset.variables.keys():
    speed = dataset.variables["surfaceSpeed"][times,:] * 3600. * 24. * 365.
elif "surfaceSpeed" not in dataset.variables.keys() and \
    all([ii in dataset.variables.keys() for ii in ['uReconstructX', 'uReconstructY']]):
        speed = np.sqrt(dataset.variables["uReconstructX"][times,:,0]**2. +
                        dataset.variables["uReconstructY"][times,:,0]**2.)
        speed *= 3600. * 24. * 365.  # convert from m/s to m/yr
else:
    print('File does not contain surfaceSpeed or uReconstructX/Y.',
          ' Skipping velocity plot.')
    plot_speed = False


if options.interp_temp:
    temperature = dataset.variables['temperature'][times,:]
    
bedTopo = dataset.variables["bedTopography"][0,:]
print('Reading bedTopography from the first time level only. If multiple',
      'times are needed, plot_transects.py will need to be updated.')
# Use coordinates from CSV file or -x -y options, but not both.
# CSV file takes precedent if both are present.
if options.coords_file is not None:
    x = []
    y = []
    with open(options.coords_file, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
     
        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))
    if [options.x_coords, options.y_coords] != [None, None]:
        print('-c and -x/-y options were both provided. Reading from ',
              f'{options.coords_file} and ignoring -x and -y settings.')
    x = np.asarray(x)
    y = np.asarray(y)
else:
    x = np.array([float(i) for i in options.x_coords.split(',')])
    y = np.array([float(i) for i in options.y_coords.split(',')])

# increase sampling to match highest mesh resolution
total_distance, = np.cumsum( np.sqrt( np.diff(x)**2. + np.diff(y)**2. ) )
n_samples = int(round(total_distance / np.min(dataset.variables["dcEdge"][:])))
x_interp = np.interp(np.linspace(0, len(x)-1, n_samples),
                     np.linspace(0, len(x)-1, len(x)), x)
y_interp = np.interp(np.linspace(0, len(y)-1, n_samples),
                     np.linspace(0, len(y)-1, len(y)), y)

d_distance = np.zeros(len(x_interp))
for ii in np.arange(1, len(x_interp)):
    d_distance[ii] = np.sqrt( (x_interp[ii] - x_interp[ii-1])**2 +
                              (y_interp[ii] - y_interp[ii-1])**2 )

distance = np.cumsum(d_distance) / 1000.  # in km for plotting

transectFig, transectAx = plt.subplots(2,1, sharex=True, layout='constrained')
thickAx = transectAx[0]
thickAx.grid()
speedAx = transectAx[1]
speedAx.grid()
timeColors = cm.plasma(np.linspace(0,1,len(times)))

plt.rcParams.update({'font.size': 16})

bed_interpolant = LinearNDInterpolator(np.vstack((xCell, yCell)).T, bedTopo)
bed_transect = bed_interpolant(np.vstack((x_interp, y_interp)).T)

for i, time in enumerate(times):
    thk_interpolant = LinearNDInterpolator(
                        np.vstack((xCell, yCell)).T, thk[i,:])
    thk_transect = thk_interpolant(np.vstack((x_interp, y_interp)).T)
    lower_surf = np.maximum( -910. / 1028. * thk_transect, bed_transect)
    lower_surf_nan = lower_surf.copy()  # for plotting
    lower_surf_nan[thk_transect==0.] = np.nan
    upper_surf = lower_surf + thk_transect
    upper_surf_nan = upper_surf.copy()  # for plotting
    upper_surf_nan[thk_transect==0.] = np.nan
    thickAx.plot(distance, lower_surf_nan, color=timeColors[i])
    thickAx.plot(distance, upper_surf_nan, color=timeColors[i])

    if plot_speed:
        speed_interpolant = LinearNDInterpolator(
                                np.vstack((xCell, yCell)).T, speed[i,:])
        speed_transect = speed_interpolant(np.vstack((x_interp, y_interp)).T)
        speed_transect[thk_transect == 0.] = np.nan
        speedAx.plot(distance, speed_transect, color=timeColors[i])

    if options.interp_temp:
        layer_thk = np.zeros((len(thk_transect), nVertLevels + 1))
        layer_midpoints = np.zeros((len(thk_transect), nVertLevels))
        layer_interfaces = np.zeros((len(thk_transect), nVertLevels + 1))
        layer_thk[:,0] = 0.
        for ii in range(len(thk_transect)):
            layer_thk[ii,1:] = np.cumsum(layerThicknessFractions *
                                         thk_transect[ii])
            layer_midpoints[i,:] = upper_surf[ii] - (layer_thk[ii,1:] +
                                                     layer_thk[ii,0:-1]) / 2.
            layer_interfaces[ii,:] = upper_surf[ii] - layer_thk[ii,:]

        temp_transect = np.zeros((len(x_interp), nVertLevels))
        for lev in range(nVertLevels):
            print(f'Interpolating temperature for level {lev}')
            temp_interpolant = LinearNDInterpolator(
                                    np.vstack((xCell, yCell)).T,
                                    temperature[i,:,lev])
            temp_transect[:, lev] = temp_interpolant(
                                        np.vstack((x_interp, y_interp)).T)

thickAx.plot(distance, bed_transect, color='black')

if options.interp_temp:
    temp_transect[temp_transect == 0.] = np.nan
    temp_plot = thickAx.pcolormesh( np.tile(distance, (nVertLevels+1,1)).T,
                                    layer_interfaces[:,:], temp_transect[1:,:],
                                    cmap='YlGnBu_r',
                                    vmin=240., vmax=273.15)

speedAx.set_xlabel('Distance (km)')
speedAx.set_ylabel('Surface\nspeed (m/yr)')
thickAx.set_ylabel('Elevation\n(m asl)')

if options.interp_temp:
    temp_cbar = plt.colorbar(temp_plot)
    temp_cbar.set_label('Temperature (K)')
    temp_cbar.ax.tick_params(labelsize=12)

if (len(times) > 1):
    time_cbar = plt.colorbar(cm.ScalarMappable(cmap='plasma'), ax=thickAx)
    time_cbar.ax.tick_params(labelsize=12)
    if use_yrs:
        time_cbar.set_label('Year')
    else:
        time_cbar.set_label('time index')
    time_cbar.set_ticks(times / np.max(times))
    time_cbar.set_ticklabels(times_list)

if options.save_filename is not None:
    transectFig.savefig(options.save_filename, dpi=400, bbox_inches='tight')

plt.show()

dataset.close()
