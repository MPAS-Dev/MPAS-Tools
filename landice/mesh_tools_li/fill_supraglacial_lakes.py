#!/usr/bin/env python3
"""
fill_supraglacial_lakes.py

fill depressions in the grounded ice-sheet upper surface
"""

import argparse
import heapq
import mosaic
import numpy as np
import sys
import xarray as xr

from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from compass.landice.mesh import mpas_flood_fill

rhoi = 910.0
rhoo = 1028.0


class CyclicNormalize(mcolors.Normalize):
    def __init__(self, cmin=0, cmax=1, vmin=0, vmax=1, clip=False):
        self.cmin = cmin
        self.cmax = cmax
        mcolors.Normalize.__init__(self, vmin, vmax, clip=clip)

    def __call__(self, value, clip=False):
        x, y = [self.cmin, self.cmax], [0, 1]
        return np.ma.masked_array(np.interp(value, x, y, period=self.cmax - self.cmin))


def parse_arguments():
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments containing input file paths, parameters, and options.
    """
    parser = argparse.ArgumentParser(
        description="Run a priority-flood fill on a 2D unstructured potential surface."
    )

    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Path to input file"
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Path to output file"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output for debugging."
    )

    return parser.parse_args()

def run(args):
    """
    Load MALI geometry data and calculates needed derived fields

    Returns
    -------
    xarray.dataset
        
        
    """
    ds = xr.open_dataset(args.input, decode_times=False, decode_cf=False)
    thk = ds.thickness[0, :]
    bed = ds.bedTopography[0, :]
    ds["grd_mask"] = ((thk * rhoi / rhoo +bed) > 0.0) * (thk > 0.0)
    ds["upperSurfaceGrd"] = ds.grd_mask * (bed + thk)
    usrf = ds['upperSurfaceGrd'].values[:]
    nCells = ds.nCells.values
    cellsOnCell = ds.cellsOnCell.values
    nEdgesOnCell = ds.nEdgesOnCell.values

    # Create mask for grounded margin, defined here to include transition
    # between floating and grounded
    ds["grd_margin_mask"] = xr.zeros_like(ds.upperSurfaceGrd, dtype=int)
    for iCell in np.argwhere(ds.grd_mask.values == 1):
        neighbor_indices = cellsOnCell[iCell, :nEdgesOnCell[iCell][0]] - 1
        neighbor_indices = neighbor_indices[neighbor_indices >= 0]
        if ds.grd_mask.values[neighbor_indices].sum() != len(neighbor_indices):
            ds.grd_margin_mask[iCell] = 1

    # create mask for hydrologically connected region to margin
    #print("Creating mask of hydrologically connected region to margin")
    #ds['mgn_draining_mask'] = xr.zeros_like(ds.upperSurfaceGrd, dtype=int)
    #iter = 0
    #seed_mask = ds['grd_margin_mask'].values[:]
    #print(f"Size of grd_margin_mask={len(seed_mask)}")
    #grow_mask = ds["grd_mask"].values[:]
    #new_keep_mask = seed_mask.copy()
    #keep_mask = new_keep_mask * 0
    #n_mask_cells = keep_mask.sum()
    #grow_iters=sys.maxsize
    #new_ind = np.nonzero((new_keep_mask - keep_mask) == 1)[0]
    #for iter in range(grow_iters):
    #    new_keep_mask = keep_mask.copy()
    #    print(f'iter={iter}, keep_mask size={keep_mask.sum()}')

    #    for iCell in new_ind:
    #        neighs = cellsOnCell[iCell, :nEdgesOnCell[iCell]] - 1
    #        neighs = neighs[neighs >= 0]  # drop garbage cell
    #        for jCell in neighs:
    #            if grow_mask[jCell] == 1 and usrf[jCell] >= usrf[iCell]:
    #                new_keep_mask[jCell] = 1

    #    # only search over new cells added previous iteration
    #    new_ind = np.nonzero((new_keep_mask - keep_mask) == 1)[0]
    #    keep_mask = new_keep_mask.copy()

    #    n_mask_cells_new = keep_mask.sum()
    #    if n_mask_cells_new == n_mask_cells:
    #        break
    #    n_mask_cells = n_mask_cells_new
    #    iter += 1
    #ds["mgn_draining_mask"][:] = keep_mask

    # ---------------
    # find depressions
    grd_mask = ds['grd_mask'].values[:]
    margin = ds['grd_margin_mask'].values[:]
    filled = usrf.copy()
    visited = np.zeros_like(usrf, dtype=bool)
    pq = []

    # add margin cells to heap to start
    for iCell in np.where(margin == 1)[0]:
        heapq.heappush(pq, (usrf[iCell], iCell))
        visited[iCell] = True

    while pq:
        ht, iCell = heapq.heappop(pq)
        neighbor_indices = cellsOnCell[iCell, 0:nEdgesOnCell[iCell]] - 1
        neighbor_indices = neighbor_indices[neighbor_indices >= 0]
        for jCell in neighbor_indices:
            if not visited[jCell] and grd_mask[jCell]:
                visited[jCell] = True
                filled[jCell] = max(usrf[jCell], ht)
                heapq.heappush(pq, (filled[jCell], jCell))
    ds['usrf_filled'] = xr.zeros_like(ds.upperSurfaceGrd)
    ds['usrf_filled'][:] = filled

    diff = ds['usrf_filled'] - ds['upperSurfaceGrd']
    diff = diff.values

    print("Filling statistics:")
    print(f"#cells adjusted by >  0m: {(diff>0.).sum()}")
    print(f"#cells adjusted by >  5m: {(diff>5.).sum()}")
    print(f"#cells adjusted by > 10m: {(diff>10.).sum()}")
    print(f"#cells adjusted by > 20m: {(diff>20.).sum()}")
    print(f"#cells adjusted by > 50m: {(diff>50.).sum()}")
    print(f"#cells adjusted by >100m: {(diff>100.).sum()}")


    dsout = xr.open_dataset(args.input, decode_times=False, decode_cf=False)
    dsout['thickness'] = dsout['thickness'] + (ds['usrf_filled'] - ds['upperSurfaceGrd'])
    dsout.to_netcdf(args.output)


    if args.verbose:
        # plot result
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(14, 10),
                                 constrained_layout=True,
                                 sharex=True, sharey=True)
        axes = axes.flatten()
        descriptor = mosaic.Descriptor(ds)

        ax = 0
        pc = mosaic.polypcolor(axes[ax], descriptor,
                               c=(thk > 0.0), aa=False)
        axes[ax].set_title(f"ice mask")

        ax += 1
        pc = mosaic.polypcolor(axes[ax], descriptor,
                               ds.grd_mask, aa=False)
        axes[ax].set_title(f"grd mask")

        ax += 1
        pc = mosaic.polypcolor(axes[ax], descriptor,
                               ds.grd_margin_mask, aa=False)
        axes[ax].set_title(f"grounded margin mask")

        #pc = mosaic.polypcolor(axes[3], descriptor,
        #                       ds.mgn_draining_mask + ds.grd_mask, aa=False)
        #axes[3].set_title(f"margin draining mask + grd mask")

        ax += 1
        period = 200.0
        norm = mcolors.Normalize(vmin=0, vmax=period)
        cyclicnorm = CyclicNormalize(cmin=0., cmax=400.0, vmin=0.0, vmax=4500.0)
        pc = mosaic.polypcolor(axes[ax], descriptor,
                               ds['upperSurfaceGrd'], aa=False,
                               cmap="hsv", norm=cyclicnorm)
        #pc = mosaic.polypcolor(axes[ax], descriptor,
        #                       ds.mgn_draining_mask * np.nan,
        #                       ec=np.where(ds.mgn_draining_mask.values == 1, "white", "black"))
        ##fig.colorbar(pc, ax=axes[ax], fraction=0.1,
        ##             label=f"surface elevation (m)")
        pc = mosaic.polypcolor(axes[ax], descriptor,
                               ds['usrf_filled'] * np.nan,
                               ec=np.where(usrf == filled, "white", "black"),
                               lw=np.where(usrf == filled, 0.01, 0.5))
        axes[ax].set_title(f"upper surface")

        ax += 1
        pc = mosaic.polypcolor(axes[ax], descriptor,
                               ds['usrf_filled'], aa=False,
                               cmap="hsv", norm=cyclicnorm)
        pc = mosaic.polypcolor(axes[ax], descriptor,
                               ds['usrf_filled'] * np.nan,
                               ec=np.where(usrf == filled, "white", "black"),
                               lw=np.where(usrf == filled, 0.01, 0.5))
        axes[ax].set_title(f"upper surface filled")

        ax += 1
        diff = ds['usrf_filled'] - ds['upperSurfaceGrd']
        diff = diff.where(diff != 0.0)
        cmap = cm.get_cmap("rainbow", 10).copy()
        cmap.set_bad(color="white")
        pc = mosaic.polypcolor(axes[ax], descriptor,
                               diff,
                               vmin=0, vmax=50,
                               aa=False, cmap=cmap)
        axes[ax].set_title(f"upper surface adjustment")
        fig.colorbar(pc, ax=axes[ax], fraction=0.1)


        for ax in axes:
            ax.set_aspect('equal')
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")




def main():
    """
    Main entry point for the script.
    """
    args = parse_arguments()

    if args.verbose:
        print("Arguments parsed successfully:")
        print(args)

    # perform processing
    run(args)

    if args.verbose:
        print("Execution complete.")
        plt.show()

    return 0


if __name__ == "__main__":
    sys.exit(main())
