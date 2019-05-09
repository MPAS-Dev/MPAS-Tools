#!/usr/bin/env python
"""
Name: widen_transect_edge_masks.py
Author: Mark Petersen, Xylar Asay-Davis

Alter transects to be at least two cells wide.  This is used for critical
passages, to avoid sea ice blockage.  Specifically, mark cells on both sides
of each transect edge mask as a water cell.
"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import argparse
import xarray

from mpas_tools.ocean.coastline_alteration import widen_transect_edge_masks


if __name__ == '__main__':

    parser = \
        argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--mask_file", dest="mask_filename",
                        help="Mask file with cell and edge transect masks.",
                        metavar="MASKFILE",
                        required=True)
    parser.add_argument("-m", "--mesh_file", dest="mesh_filename",
                        help="MPAS Mesh filename.", metavar="MESHFILE",
                        required=True)
    parser.add_argument("-o", "--out_file", dest="out_filename",
                        help="Output mask file,different from input filename.",
                        metavar="MASKFILE",
                        required=True)
    parser.add_argument("-l", "--latitude_threshold",
                        dest="latitude_threshold",
                        help="Minimum latitude, degrees, for transect "
                             "widening.",
                        required=False, type=float, default=43.0)
    args = parser.parse_args()

    dsMask = xarray.open_dataset(args.mask_filename)

    dsMesh = xarray.open_dataset(args.mesh_filename)

    dsMask = widen_transect_edge_masks(dsMask, dsMesh, args.latitude_threshold)
    dsMask.to_netcdf(args.out_filename)
