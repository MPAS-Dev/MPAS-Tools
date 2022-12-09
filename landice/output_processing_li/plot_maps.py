#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 9, 2022

@author: Trevor Hillebrand, Matthew Hoffman
"""
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-r", dest="runs", help="path to .nc file or dir containing output.nc file (strings separated by commas; no spaces)", default=None, metavar="FILENAME")
parser.add_option("-t", dest="timeLevels", help="integer time levels at which to plot (int separated by commas; no spaces)", default='-1')
parser.add_option("-v", dest="variables", help="variable(s) to plot (list separated by commas; no spaces)", default='thickness')
parser.add_option("-c", dest="colormaps", help="colormaps to use for plotting (list separated by commas, no spaces", default=None)
parser.add_option("-s", dest="saveNames", help="filename for saving. If empty or None, will plot to screen instead of saving.", default=None, metavar="FILENAME")

options, args = parser.parse_args()
runs = options.runs.split(',') # split run directories into list
variables = options.variables.split(',')
timeLevs = options.timeLevels.split(',')  # split time levels into list

if options.colormaps is not None:
    colormaps = options.colormaps.split(',')
else:
    colormaps = ['viridis'] * len(variables)

if options.saveNames is not None:
    saveNames = options.saveNames.split(',')

initialExtentValue = 1
dynamicValue = 2
floatValue = 4
groundingLineValue = 256

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
    # If would be more efficient to do this outside
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
            axs.append(plt.subplot(gs[run][row, col]))

    varPlot[run] = {}  # is a dict of dicts too complicated?
    cbars = []
    # Loop over variables
    for row, (variable, colormap, cbar_ax) in enumerate(zip(variables, colormaps, cbar_axs)):
        var_to_plot = f.variables[variable][:]
        units = f.variables[variable].units
        varPlot[run][variable] = []
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
            timeLev = int(timeLev) #these are strings for some reason; make int to index
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
                                  triang, var_to_plot[timeLev,:], cmap=colormap,
                                  shading='flat', vmin=np.nanquantile(var_to_plot, 0.01),
                                  vmax=np.nanquantile(var_to_plot, 0.99)))
            axs[index].set_aspect('equal')
            axs[index].set_title(f'year = {yr[timeLev]:0.2f}')

        cbars.append(Colorbar(ax=cbar_ax, mappable=varPlot[run][variable][0], orientation='vertical',
                 label=f'{variable} (${units}$)'))

    if options.saveNames is not None:
        figs[run].savefig(saveNames[ii], dpi=400, bbox_inches='tight')
    
    f.close()

plt.show()
