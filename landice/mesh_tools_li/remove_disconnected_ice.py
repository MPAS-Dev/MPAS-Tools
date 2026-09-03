#!/usr/bin/env python
'''
Remove small, dynamically disconnected clusters of ice from a MALI mesh.

Ice cells (thickness greater than a threshold) are grouped into connected
components using the mesh cell adjacency (cellsOnCell).  A flood-fill is
seeded at a user-defined location (default: the centre of the mesh) to
identify the main ice body; the component containing the seed cell is always
preserved.  Every other ice cluster that is not connected to the main body is
considered disconnected.  Any disconnected cluster whose size is below a
user-defined threshold (default: 10 cells) has its ice thickness set to zero.
In addition, disconnected clusters that are entirely floating are removed
regardless of size (detached icebergs / ice-shelf fragments).  Larger grounded
disconnected clusters (e.g. genuine separate ice caps) are left unchanged.

The modified mesh is written to the required output file (--out_file), which
must differ from the input; the input file is left unchanged.

Trevor Hillebrand, 2026
'''

import os
import subprocess
import numpy as np
import xarray as xr
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from argparse import ArgumentParser
from datetime import datetime
import sys

rho_ice = 910.0
rho_sw = 1028.0


def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--file", dest="mesh_file", required=True,
                        metavar="FILENAME",
                        help="MALI mesh file to process")
    parser.add_argument("-o", "--out_file", dest="out_file", required=True,
                        metavar="FILENAME",
                        help="Output file to write (must differ from the input "
                             "file). The input is read, modified in memory, and "
                             "written here; the input file is left unchanged.")
    parser.add_argument("-x", "--x_seed", dest="x_seed", type=float,
                        default=None,
                        help="x coordinate [m] of the flood-fill seed "
                             "(default: centre of the mesh)")
    parser.add_argument("-y", "--y_seed", dest="y_seed", type=float,
                        default=None,
                        help="y coordinate [m] of the flood-fill seed "
                             "(default: centre of the mesh)")
    parser.add_argument("-n", "--min_cells", dest="min_cells", type=int,
                        default=10,
                        help="Disconnected ice clusters smaller than this many "
                             "cells are removed (default: 10)")
    parser.add_argument("-t", "--thickness_threshold",
                        dest="thickness_threshold", type=float, default=0.0,
                        help="Ice is defined as thickness greater than this "
                             "value [m] (default: 0.0)")
    parser.add_argument("--time_level", dest="time_level", type=int, default=0,
                        help="Time index of the thickness field to modify "
                             "(default: 0)")
    return parser.parse_args()


def ice_connected_components(ice_mask, cells_on_cell, n_edges_on_cell):
    '''Label connected components of the ice mask over the mesh graph.

    Returns (labels, sizes) where labels[i_cell] is the component id of each
    cell and sizes[label] is the number of cells in that component. Only
    ice-ice adjacencies create edges, so every ice cluster is its own
    component and non-ice cells are isolated singletons.'''
    n_cells = len(ice_mask)
    max_edges = cells_on_cell.shape[1]
    src = np.repeat(np.arange(n_cells), max_edges)
    dst = (cells_on_cell - 1).ravel()
    # Drop phantom neighbours (0 in cellsOnCell -> -1) and edges beyond the
    # actual edge count of each cell.
    within_edges = (np.arange(max_edges)[None, :] <
                    n_edges_on_cell[:, None]).ravel()
    valid = (dst >= 0) & within_edges
    src, dst = src[valid], dst[valid]
    # Keep only edges connecting two ice cells.
    ice_edge = ice_mask[src] & ice_mask[dst]
    src, dst = src[ice_edge], dst[ice_edge]
    graph = coo_matrix((np.ones(src.size, dtype=np.int8), (src, dst)),
                       shape=(n_cells, n_cells))
    _, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels)
    return labels, sizes


def main():
    args = parse_args()

    assert os.path.abspath(args.out_file) != os.path.abspath(args.mesh_file), \
        "Output file must be different from the input file."

    ds = xr.open_dataset(args.mesh_file, decode_times=False, decode_cf=False)
    x_cell = ds['xCell'].values
    y_cell = ds['yCell'].values
    cells_on_cell = ds['cellsOnCell'].values
    n_edges_on_cell = ds['nEdgesOnCell'].values
    thickness = ds['thickness'].isel(Time=args.time_level).values.copy()

    ice_mask = thickness > args.thickness_threshold
    n_ice = int(ice_mask.sum())
    if n_ice == 0:
        ds.close()
        sys.exit("ERROR: no ice cells (thickness > "
                 f"{args.thickness_threshold}) found in the mesh.")

    # Grounded ice satisfies the flotation criterion; a disconnected cluster
    # with no grounded cells is floating and is removed regardless of size.
    if 'bedTopography' in ds.variables:
        bed = ds['bedTopography'].isel(Time=args.time_level).values
        grounded_mask = ice_mask & (thickness * rho_ice / rho_sw + bed > 0.0)
    else:
        grounded_mask = None
        print("WARNING: 'bedTopography' not found; cannot identify floating "
              "ice, so only the size threshold will be applied.")

    # Locate the seed cell (nearest cell to the requested location; default is
    # the centre of the mesh bounding box).
    x_seed = args.x_seed if args.x_seed is not None else \
        0.5 * (x_cell.min() + x_cell.max())
    y_seed = args.y_seed if args.y_seed is not None else \
        0.5 * (y_cell.min() + y_cell.max())
    seed_cell = int(np.argmin((x_cell - x_seed) ** 2 + (y_cell - y_seed) ** 2))
    if not ice_mask[seed_cell]:
        # Snap to the nearest ice cell so the main body is well defined.
        ice_idx = np.where(ice_mask)[0]
        nearest = ice_idx[np.argmin((x_cell[ice_idx] - x_seed) ** 2 +
                                    (y_cell[ice_idx] - y_seed) ** 2)]
        print(f"WARNING: seed location ({x_seed:.1f}, {y_seed:.1f}) is on an "
              f"ice-free cell; snapping to nearest ice cell {nearest}.")
        seed_cell = int(nearest)

    labels, sizes = ice_connected_components(ice_mask, cells_on_cell,
                                             n_edges_on_cell)
    main_label = labels[seed_cell]
    print(f"Seed cell {seed_cell} at ({x_cell[seed_cell]:.1f}, "
          f"{y_cell[seed_cell]:.1f}); main ice body has "
          f"{sizes[main_label]} cells")

    # Flag clusters that contain any grounded ice.
    cluster_has_grounded = np.zeros(sizes.size, dtype=bool)
    if grounded_mask is not None:
        cluster_has_grounded[labels[grounded_mask]] = True
    else:
        cluster_has_grounded[:] = True

    # Remove disconnected clusters that are small OR entirely floating.
    ice_labels = np.unique(labels[ice_mask])
    disconnected = ice_labels[ice_labels != main_label]
    is_small = sizes[disconnected] < args.min_cells
    is_floating = ~cluster_has_grounded[disconnected]
    remove_labels = disconnected[is_small | is_floating]
    print(f"Found {len(disconnected)} disconnected ice cluster(s); removing "
          f"{len(remove_labels)} (floating: {int(is_floating.sum())}, "
          f"below {args.min_cells}-cell threshold: {int(is_small.sum())})")

    remove_mask = ice_mask & np.isin(labels, remove_labels)
    n_removed = int(remove_mask.sum())
    thickness[remove_mask] = 0.0
    ds['thickness'][args.time_level, :] = thickness
    print(f"Zeroed ice thickness on {n_removed} cell(s) in "
          f"{len(remove_labels)} disconnected cluster(s)")

    # Update the history attribute.
    this_command = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + \
        " ".join(sys.argv[:])
    if 'history' in ds.attrs:
        ds.attrs['history'] = "\n".join([this_command, ds.attrs['history']])
    else:
        ds.attrs['history'] = this_command

    # Writing NetCDF-3 directly from xarray is slow, so write NETCDF4 first
    # and convert to NETCDF3_64BIT with ncks, which is much faster.
    if 'Time' in ds.dims:
        ds.encoding['unlimited_dims'] = {'Time'}
    tmp_file = f"{args.out_file}.netcdf4.tmp"
    ds.to_netcdf(tmp_file, format='NETCDF4')
    ds.close()
    subprocess.check_call(['ncks', '-O', '-6', tmp_file, args.out_file])
    os.remove(tmp_file)
    print(f"Wrote '{args.out_file}'")


if __name__ == '__main__':
    main()
