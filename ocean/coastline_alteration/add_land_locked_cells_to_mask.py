#!/usr/bin/env python
"""
Name: add_land_locked_cells_to_mask.py
Author: Mark Petersen, Adrian Turner, Xylar Asay-Davis

Find ocean cells that are land-locked, and alter the cell
mask so that they are counted as land cells.
"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import argparse
import xarray

from mpas_tools.ocean.coastline_alteration import add_land_locked_cells_to_mask

if __name__ == '__main__':
    parser = \
        argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--input_mask_file", dest="input_mask_filename",
                        help="Mask file that includes cell and edge masks.",
                        metavar="INPUTMASKFILE", required=True)
    parser.add_argument("-o", "--output_mask_file",
                        dest="output_mask_filename",
                        help="Mask file that includes cell and edge masks.",
                        metavar="OUTPUTMASKFILE", required=True)
    parser.add_argument("-m", "--mesh_file", dest="mesh_filename",
                        help="MPAS Mesh filename.", metavar="MESHFILE",
                        required=True)
    parser.add_argument("-l", "--latitude_threshold",
                        dest="latitude_threshold",
                        help="Minimum latitude, in degrees, for transect "
                             "widening.",
                        required=False, type=float, default=43.0)
    parser.add_argument("-n", "--number_sweeps", dest="nSweeps",
                        help="Maximum number of sweeps to search for "
                             "land-locked cells.",
                        required=False, type=int, default=10)
    args = parser.parse_args()

    dsMask = xarray.open_dataset(args.input_mask_filename)

    dsMesh = xarray.open_dataset(args.mesh_filename)

    dsMask = add_land_locked_cells_to_mask(dsMask, dsMesh,
                                           args.latitude_threshold,
                                           args.nSweeps)
    dsMask.to_netcdf(args.output_mask_filename)
