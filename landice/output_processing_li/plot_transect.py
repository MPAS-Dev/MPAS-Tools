#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot bed and surface topography, surface speed, and
(optionally) temperature and thermal forcing from MPAS netCDF along a transect.
@author: trevorhillebrand, edited by cashafer
"""

import numpy as np
import csv
from netCDF4 import Dataset
from optparse import OptionParser
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
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

parser.add_option("--tf", dest="interp_thermal_forcing",
                  action="store_true", help="interpolate thermal forcing")

parser.add_option("--nofill", dest="nofill",
                  action="store_true", help="disable fill for glacier and topography \
                                             (useful for debugging TF extrapolation)")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()


# Get list of times for plotting
times_list = [i for i in options.times.split(',')]  # list of string times for plotting
times = [int(i) for i in options.times.split(',')]  # list of integer time indices

# Load MPAS netCDF file
dataset = Dataset(options.data_file, 'r')
dataset.set_always_mask(False)

# Get x and y coordinates
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

# Determine bounding box using max and min of user-provided coordinates
xCell = dataset.variables["xCell"][:]
yCell = dataset.variables["yCell"][:]
dcEdge = dataset.variables["dcEdge"][:]
x_max = np.max(x)
x_min = np.min(x)
y_max = np.max(y)
y_min = np.min(y)
indices1 = np.where(np.logical_and(
                    np.logical_and(xCell>=x_min, xCell<=x_max),
                    np.logical_and(yCell>=y_min, yCell<=y_max)) )[0]
pad = dcEdge[indices1].max() * 1.5  # pad based on the largest cell spacing in the region of interest
print(f'Using padding around region of interest of {pad} m')

x_max += pad
x_min -= pad
y_max += pad
y_min -= pad
indices = np.where( np.logical_and(
                    np.logical_and(xCell>=x_min, xCell<=x_max),
                    np.logical_and(yCell>=y_min, yCell<=y_max)) )[0]
xCell = xCell[indices]
yCell = yCell[indices]

tri_mesh = Delaunay( np.vstack((xCell, yCell)).T )

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

# Cannot plot temperature or thermal forcing for more than one time index.
if (options.interp_temp or options.interp_thermal_forcing) and (len(times) > 1):
    print('Cannot plot temperature or thermal forcing for more than one time index.' +
          ' Skipping temperature or TF interpolation and plotting.')
    options.interp_temp = False
    options.interp_thermal_forcing = False

li_mask_ValueDynamicIce = 2
cellMask = dataset.variables['cellMask'][times, indices]
cellMask_dynamicIce = (cellMask & li_mask_ValueDynamicIce) // li_mask_ValueDynamicIce
# only take thickness of dynamic ice
thk = dataset.variables["thickness"][times, indices] * cellMask_dynamicIce

plot_speed = True
# Include speed on non-dynamic ice to avoid interpolation artifacts.
if "surfaceSpeed" in dataset.variables.keys():
    speed = dataset.variables["surfaceSpeed"][times, indices] * 3600. * 24. * 365.
elif "surfaceSpeed" not in dataset.variables.keys() and \
    all([ii in dataset.variables.keys() for ii in ['uReconstructX', 'uReconstructY']]):
        speed = np.sqrt(dataset.variables["uReconstructX"][times, indices, 0]**2. +
                        dataset.variables["uReconstructY"][times, indices, 0]**2.)
        speed *= 3600. * 24. * 365.  # convert from m/s to m/yr
else:
    print('File does not contain surfaceSpeed or uReconstructX/Y.',
          ' Skipping velocity plot.')
    plot_speed = False

if options.interp_temp:
    temperature = dataset.variables['temperature'][times, indices]

if options.interp_thermal_forcing:
    # Extrapolated Ocean thermal forcing
    thermal_forcing = dataset.variables['ismip6shelfMelt_3dThermalForcing'][times, indices]
    thermal_forcing[thermal_forcing > 1.0e3] = np.nan # Large invalid TF values set to NaN
    # Number of ocean layers
    nISMIP6OceanLayers = dataset.dimensions['nISMIP6OceanLayers'].size  
    # Depths associated with thermal forcing field
    if 'ismip6shelfMelt_zOcean' in dataset.variables:
        ismip6shelfMelt_zOcean = dataset.variables['ismip6shelfMelt_zOcean'][:]
    else:
        ismip6shelfMelt_zOcean = np.linspace(-30.0, -1770.0, num=30)
        print('WARNING: ismip6shelfMelt_zOcean not found in file. Assuming values:')
        print(ismip6shelfMelt_zOcean)
    
bedTopo = dataset.variables["bedTopography"][0, indices]
print('Reading bedTopography from the first time level only. If multiple',
      'times are needed, plot_transects.py will need to be updated.')


# increase sampling to match highest mesh resolution
total_distance = np.cumsum( np.sqrt( np.diff(x)**2. + np.diff(y)**2. ) )
n_samples = int(round(total_distance[-1] / np.min(dcEdge[indices])))
x_interp = np.interp(np.linspace(0, len(x)-1, n_samples),
                     np.linspace(0, len(x)-1, len(x)), x)
y_interp = np.interp(np.linspace(0, len(y)-1, n_samples),
                     np.linspace(0, len(y)-1, len(y)), y)

d_distance = np.zeros(len(x_interp))
for ii in np.arange(1, len(x_interp)):
    d_distance[ii] = np.sqrt( (x_interp[ii] - x_interp[ii-1])**2 +
                              (y_interp[ii] - y_interp[ii-1])**2 )

distance = np.cumsum(d_distance) / 1000.  # in km for plotting

transectFig, transectAx = plt.subplots(2,1, sharex=True, layout='constrained', figsize=(8,6))
thickAx = transectAx[0]
thickAx.grid()
speedAx = transectAx[1]
speedAx.grid()
timeColors = cm.plasma(np.linspace(0,1,len(times)))

bed_interpolant = LinearNDInterpolator(tri_mesh, bedTopo)
bed_transect = bed_interpolant(np.vstack((x_interp, y_interp)).T)

for i, time in enumerate(times):
    thk_interpolant = LinearNDInterpolator(tri_mesh, thk[i,:])
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
        speed_interpolant = LinearNDInterpolator(tri_mesh, speed[i,:])
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
            temp_interpolant = LinearNDInterpolator(tri_mesh, temperature[i,:,lev])
            temp_transect[:, lev] = temp_interpolant(
                                        np.vstack((x_interp, y_interp)).T)
 
    if options.interp_thermal_forcing:
        ocean_layer_interfaces = np.zeros((len(thk_transect), nISMIP6OceanLayers + 1))  # Ocean layer depths
        ocean_layer_interfaces[:,0] = 0.        # First interface is the topmost z = 0 sea level interface

        # For the entire length of the transect, set the layer interface depths
        for ii in range(len(thk_transect)):
            ocean_layer_interfaces[ii,1:] = ismip6shelfMelt_zOcean[:] + ismip6shelfMelt_zOcean[0]

        # Create the thermal forcing transect using the TF interpolator at each ocean lev
        thermal_forcing_transect = np.zeros((len(x_interp), nISMIP6OceanLayers))

        for ocean_lev in range(nISMIP6OceanLayers):
            print(f'Interpolating ocean thermal forcing for ocean level {ocean_lev}')
            thermal_forcing_interpolant = LinearNDInterpolator(tri_mesh, thermal_forcing[i, :, ocean_lev])
            thermal_forcing_transect[:, ocean_lev] = thermal_forcing_interpolant(np.vstack( (x_interp, y_interp) ).T)

thickAx.plot(distance, bed_transect, color='black')

if options.interp_temp:
    temp_transect[temp_transect == 0.] = np.nan
    temp_plot = thickAx.pcolormesh( np.tile(distance, (nVertLevels+1,1)).T,
                                    layer_interfaces[:,:], temp_transect[1:,:],
                                    cmap='YlGnBu_r',
                                    vmin=240., vmax=273.15)

if options.interp_thermal_forcing:
    thermal_forcing_plot = thickAx.pcolormesh(np.tile(distance, (nISMIP6OceanLayers+1, 1)).T,
                                               ocean_layer_interfaces[:, :], thermal_forcing_transect[1:, :],
                                               cmap='YlOrRd', vmin=0, vmax=5.5)
    thickAx.grid(False)
    # to avoid always plotting to the deepest ocean z-level,
    # manually set a reasonable range
    thickAx.set_ylim([bedTopo.min() - 50.0, upper_surf.max() + 50.0])

# add some fill
if not options.nofill:
    if not options.interp_temp and len(times) == 1:
        # only color glacier if not plotting temperature and if there is only one time
        thickAx.fill_between(distance, upper_surf_nan, lower_surf_nan, color='xkcd:ice blue')
    # always fill beneath bed
    thickAx.fill_between(distance, bed_transect, thickAx.get_ylim()[0], color='xkcd:greyish brown')

thickAx.tick_params(axis='x', labelbottom=True)
speedAx.set_xlabel('Distance (km)')
thickAx.set_xlabel('Distance (km)')
speedAx.set_ylabel('Surface\nspeed (m/yr)')
thickAx.set_ylabel('Elevation\n(m asl)')

if options.interp_temp:
    temp_cbar = plt.colorbar(temp_plot)
    temp_cbar.set_label('Temperature (K)')

if options.interp_thermal_forcing:
    tf_cbar = plt.colorbar(thermal_forcing_plot)
    tf_cbar.set_label('Thermal Forcing (C)')

if (len(times) > 1):
    time_cbar = plt.colorbar(cm.ScalarMappable(cmap='plasma'), ax=thickAx)
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

