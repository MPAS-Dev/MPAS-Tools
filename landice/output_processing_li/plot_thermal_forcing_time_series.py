#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot time series of thermal forcing at specified
locations and depths. 
@author: trevorhillebrand
"""


import numpy as np
import csv
from netCDF4 import Dataset
from optparse import OptionParser
from scipy.interpolate import LinearNDInterpolator, interp1d
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


parser = OptionParser(description='Plot transect from MPAS netCDF')
parser.add_option("-f", "--file", dest="thermal_forcing_file",
                  help="List of MPAS netCDF files that contains the ismip6shelfMelt_3dThermalForcing" \
                       " field and zOcean variable. Comma-separated, no spaces.")
parser.add_option("-m", "--mesh", dest="mesh_file",
                  help="the MPAS netCDF file that contains the mesh variable, as well as thickness and bedTopography")
parser.add_option("--start_time", dest="start_time", default="0",
                  help="beginning of time range to plot")
parser.add_option("--end_time", dest="end_time", default="-1",
                  help="end of time range to plot")
parser.add_option('-c', dest='coords_file', default=None,
                  help='CSV file containing x in first column, y in second. No header.')
parser.add_option('-x', dest='x_coords', default=None,
                  help='List of x coordinates of transect if not \
                        using a csv file. Comma-separated, no spaces.')
parser.add_option('-y', dest='y_coords', default=None,
                  help='List of y coordinates of transect if not \
                        using a csv file. Comma-separated, no spaces.')
parser.add_option('-d', dest='depth', default=None,
                  help='Depth in meters at which to plot thermal forcing.' \
                        ' If a single value, the script will use linear 1d' \
                        ' interpolation to determine the thermal forcing' \
                        ' at that depth. If two values are given, the script' \
                        ' will provide the average over that depth range.')
parser.add_option("--seafloor", dest="plot_seafloor_thermal_forcing",
                  action="store_true",
                  help="plot thermal forcing at the seafloor, instead of at specific depth")
parser.add_option("--shelf_base", dest="plot_shelf_base_thermal_forcing",
                  action="store_true",
                  help="plot thermal forcing at the base of floating ice, instead of at specific depth")
parser.add_option("--average", dest="plot_average_thermal_forcing",
                  action="store_true", help='Whether to plot average' \
                  ' thermal forcing across all coordinates provided.')
parser.add_option("--n_samples", dest="n_samples", default=None,
                  help="Number of random samples to take from provided coordinates.")
parser.add_option("-s", dest="save_filename", default=None,
                   help="File to save figure to.")

options, args = parser.parse_args()

rhoi = 910.
rhosw = 1028.
start_year = 1995  # first year in TF forcings

times_list = [options.start_time, options.end_time]  # list of string times for plotting
times = [int(i) for i in times_list]  # list of integer time indices

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

if options.n_samples is not None:
   rand_idx = np.random.choice(x.shape[0], int(options.n_samples), replace=False)
   print(f"Using {options.n_samples} random samples from the {x.shape[0]} points provided.")
   x = x[rand_idx]
   y = y[rand_idx]

# Mesh and geometry fields
mesh = Dataset(options.mesh_file, 'r')
mesh.set_always_mask(False)
bed = mesh.variables["bedTopography"][:]
thk = mesh.variables["thickness"][:]
xCell = mesh.variables["xCell"][:]
yCell = mesh.variables["yCell"][:]
nCells = mesh.dimensions['nCells'].size
areaCell = mesh.variables["areaCell"][:]
mesh.close()

fig, ax = plt.subplots(1,1, layout='constrained')

def interp_and_plot(tf_file, plot_ax):
    # Thermal forcing fields
    tf_data = Dataset(tf_file, 'r')
    tf_data.set_always_mask(False)
    tf = tf_data.variables["ismip6shelfMelt_3dThermalForcing"][:]
    z = tf_data.variables["ismip6shelfMelt_zOcean"][:]
    n_time_levs = tf_data.dimensions["Time"].size
    tf_data.close()
    if times[1] == -1:
        times[1] = n_time_levs - 1
    plot_times = np.arange(times[0], times[1], step=1)  # annual posting
    
    # Find nearest cell to desired x,y locations
    nearest_cell = []
    for x_i, y_i in zip(x, y):
        nearest_cell.append(np.argmin( np.sqrt((xCell - x_i)**2. + (yCell - y_i)**2) ))
    
    # Find depth to seafloor or ice-shelf base
    if options.depth is not None:
        depth = [float(i) for i in options.depth]
    if options.plot_seafloor_thermal_forcing:
        depth = bed[0, nearest_cell]
    elif options.plot_shelf_base_thermal_forcing:
        # Assume this is floating ice
        depth = -1. * rhoi / rhosw * thk[0, nearest_cell]

    # Clip depth to within the range of the TF data
    depth[depth > np.max(z)] = np.max(z)
    depth[depth < np.min(z)] = np.min(z)
    
    # Vertical interpolation of ocean forcing.
    tf_depth = []
    for time in plot_times:
        tf_tmp = []
        for cell, cell_depth in zip(nearest_cell, depth):
            tf_interp = interp1d(z, tf[time, cell, :])
            tf_tmp.append(tf_interp(cell_depth))
        tf_depth.append(tf_tmp)

    if "UKESM" in tf_file:
        plot_color = 'tab:blue'
    else:
        plot_color = 'tab:grey'
 
    if options.plot_average_thermal_forcing:
        tf_avg = np.mean(tf_depth, axis=1)
        tf_std = np.std(tf_depth, axis=1)
        plot_ax.plot(plot_times + start_year, tf_avg, c=plot_color)
        plot_ax.fill_between(plot_times + start_year, tf_avg - tf_std,
                             tf_avg + tf_std, fc=plot_color,
                             alpha = 0.5)
    else:
        plot_ax.plot(plot_times + start_year, tf_depth, c=plot_color)

tf_files = [i for i in options.thermal_forcing_file.split(',')]
for tf_file in tf_files:
    print(f"Plotting from {tf_file}")
    interp_and_plot(tf_file, ax)
ax.set_xlabel("Year")
ax.set_ylabel("Thermal Forcing (°C)")
ax.grid('on')
if options.save_filename is not None:
    fig.savefig(options.save_filename, dpi=400, bbox_inches='tight')    
plt.show()
