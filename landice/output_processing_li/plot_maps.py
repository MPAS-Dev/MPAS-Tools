#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 9, 2022

@author: Trevor Hillebrand, Matthew Hoffman

Script to plot snapshot maps of MALI output for an arbitrary number of files,
variables, and output times. There is no requirement for all output files
to be on the same mesh. Each output file gets its own figure, each
variable gets its own row, and each time gets its own column. Three contours
are automatically plotted, showing intial ice extent (black), initial
grounding-line position (grey), and grounding-line position at the desired
time (white).

"""
import numpy as np
from netCDF4 import Dataset
import argparse
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
from matplotlib.colors import Normalize, TwoSlopeNorm


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-r", dest="runs", default=None, metavar="FILENAME",
                    help="path to .nc file or dir containing output.nc \
                          file (strings separated by commas; no spaces)")
parser.add_argument("-t", dest="timeLevels", default="-1",
                    help="integer time levels at which to plot \
                          (int separated by commas; no spaces)")
parser.add_argument("-v", dest="variables", default='thickness',
                    help="variable(s) to plot (list separated by commas; no spaces)")
parser.add_argument("-l", dest="log_plot", default=None,
                    help="Whether to plot the log10 of each variable \
                          (True or False list separated by commas; no spaces)")
parser.add_argument("-c", dest="colormaps", default=None,
                    help="colormaps to use for plotting (list separated by commas \
                          , no spaces). This overrides default colormaps.")
parser.add_argument("-s", dest="saveNames", default=None, metavar="FILENAME",
                    help="filename for saving. If empty or None, will plot \
                          to screen instead of saving.")

args = parser.parse_args()
runs = args.runs.split(',') # split run directories into list
variables = args.variables.split(',')
timeLevs = args.timeLevels.split(',')  # split time levels into list
# convert timeLevs to list of ints
timeLevs = [int(i) for i in timeLevs]

if args.log_plot is not None:
    log_plot = args.log_plot.split(',')
else:
    log_plot = [False] * len(variables)

if args.colormaps is not None:
    colormaps = args.colormaps.split(',')
else:
    colormaps = ['viridis'] * len(variables)

if args.saveNames is not None:
    saveNames = args.saveNames.split(',')

# Set up a dictionary of default colormaps for common variables.
# These can be overridden by the -c flag.
defaultColors = {'thickness' : 'Blues',
                 'surfaceSpeed' : 'plasma',
                 'basalSpeed' : 'plasma',
                 'bedTopography' : 'BrBG',
                 'floatingBasalMassBalApplied' : 'cividis'
                }

if args.colormaps is not None:
    colormaps = args.colormaps.split(',')
else:
    colormaps = []
    for variable in variables:
        if variable in defaultColors.keys():
            colormaps.append(defaultColors[variable])
        else:
            # All other variables default to viridis
            colormaps.append('viridis')

# Set bitmask values
initialExtentValue = 1
dynamicValue = 2
floatValue = 4
groundingLineValue = 256

# List of diverging colormaps for use in plotting bedTopography.
# I don't see a way around hard-coding this.
divColorMaps = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

def dist(i1, i2, xCell, yCell):  # helper distance fn
    dist = ((xCell[i1]-xCell[i2])**2 + (yCell[i1]-yCell[i2])**2)**0.5
    return dist

# Loop over runs
# Each run gets its own figure
# Each variable gets its own row
# Each time level gets its own column
varPlot = {}
figs = {}
gs = {}
for ii, run in enumerate(runs):
    if '.nc' not in run:
        run = run + '/output.nc'
    f = Dataset(run, 'r')
    yr = f.variables['daysSinceStart'][:] / 365.
    f.set_auto_mask(False)

    # Get mesh geometry and calculate triangulation. 
    # It would be more efficient to do this outside
    # this loop if all runs are on the same mesh, but we
    # want this to be as general as possible.
    xCell = f.variables["xCell"][:]
    yCell = f.variables["yCell"][:]
    dcEdge = f.variables["dcEdge"][:]

    triang = tri.Triangulation(xCell, yCell)
    triMask = np.zeros(len(triang.triangles))
    maxDist = np.max(dcEdge) * 2.0  # maximum distance in m of edges between points. Make twice dcEdge to be safe
    for t in range(len(triang.triangles)):
        thisTri = triang.triangles[t, :]
        if dist(thisTri[0], thisTri[1], xCell, yCell) > maxDist:
            triMask[t] = True
        if dist(thisTri[1], thisTri[2], xCell, yCell) > maxDist:
            triMask[t] = True
        if dist(thisTri[0], thisTri[2], xCell, yCell) > maxDist:
            triMask[t] = True
    triang.set_mask(triMask)

    # set up figure for this run
    figs[run] = plt.figure(figsize=(15,7))
    figs[run].suptitle(run)
    nRows = len(variables)
    nCols = len(timeLevs) + 1

    # last column is for colorbars
    gs[run] = gridspec.GridSpec(nRows, nCols,
                           height_ratios=[1] * nRows,
                           width_ratios=[1] * (nCols - 1) + [0.1])
    axs = []
    cbar_axs = []
    for row in np.arange(0, nRows):
        cbar_axs.append(plt.subplot(gs[run][row,-1]))
        for col in np.arange(0, nCols-1):
            if axs == []:
                axs.append(plt.subplot(gs[run][row, col]))
            else:
                axs.append(plt.subplot(gs[run][row, col], sharex=axs[0], sharey=axs[0]))

    varPlot[run] = {}  # is a dict of dicts too complicated?
    cbars = []
    # Loop over variables
    for row, (variable, log, colormap, cbar_ax) in enumerate(zip(variables, log_plot, colormaps, cbar_axs)):
        var_to_plot = f.variables[variable][:]

        if log == 'True':
            var_to_plot = np.log10(var_to_plot)
            # Get rid of +/- inf values that ruin vmin and vmax
            # calculations below.
            var_to_plot[np.isinf(var_to_plot)] = np.nan
            colorbar_label_prefix = 'log10 '
        else:
            colorbar_label_prefix = ''
        units = f.variables[variable].units
        varPlot[run][variable] = []

        # Plot bedTopography on an asymmetric colorbar if appropriate
        if ( (variable == 'bedTopography') and
             (np.nanquantile(var_to_plot[timeLevs, :], 0.99) > 0.) and
             (colormap in divColorMaps) ):
           norm = TwoSlopeNorm(vmin=np.nanquantile(var_to_plot[timeLevs, :], 0.01),
                               vmax=np.nanquantile(var_to_plot[timeLevs, :], 0.99),
                               vcenter=0.)
        else:
            norm = Normalize(vmin=np.nanquantile(var_to_plot[timeLevs, :], 0.01),
                             vmax=np.nanquantile(var_to_plot[timeLevs, :], 0.99))

        if 'cellMask' in f.variables.keys():
            cellMask = f.variables["cellMask"][:]
            floatMask = (cellMask & floatValue) // floatValue
            dynamicMask = (cellMask & dynamicValue) // dynamicValue
            groundingLineMask = (cellMask & groundingLineValue) // groundingLineValue
            initialExtentMask = (cellMask & initialExtentValue) // initialExtentValue
        else:
            print(f'cellMask is not present in output file {run}')
        
        # Loop over time levels
        for col, timeLev in enumerate(timeLevs):
            index = row * (nCols - 1) + col
            # plot initial grounding line position, initial extent, and GL position at t=timeLev
            axs[index].tricontour(triang, groundingLineMask[0, :],
                              levels=[0.9999], colors='grey', linestyles='solid')
            axs[index].tricontour(triang, groundingLineMask[timeLev, :],
                              levels=[0.9999], colors='white', linestyles='solid')
            axs[index].tricontour(triang, initialExtentMask[timeLev, :],
                              levels=[0.9999], colors='black', linestyles='solid')

            # Plot 2D field at each desired time. Use quantile range of 0.01-0.99 to cut out
            # outliers. Could improve on this by accounting for areaCell, as currently all cells
            # are weighted equally in determining vmin and vmax.
            varPlot[run][variable].append(
                              axs[index].tripcolor(
                                  triang, var_to_plot[timeLev, :], cmap=colormap,
                                  shading='flat', norm=norm))
            axs[index].set_aspect('equal')
            axs[index].set_title(f'year = {yr[timeLev]:0.2f}')

        cbars.append(Colorbar(ax=cbar_ax, mappable=varPlot[run][variable][0], orientation='vertical',
                 label=f'{colorbar_label_prefix}{variable} (${units}$)'))

    if args.saveNames is not None:
        figs[run].savefig(saveNames[ii], dpi=400, bbox_inches='tight')
    
    f.close()

plt.show()
