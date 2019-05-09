#!/usr/bin/env python
"""
Name: add_critical_land_blockages_to_mask.py
Author: Xylar Asay-Davis

Add transects that identify critical regions where narrow strips of land block
ocean flow.  These are, essentially, the opposite of critical passages, which
must remain open for ocean flow.
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import argparse

from mpas_tools.ocean.coastline_alteration import add_critical_land_blockages


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
    parser.add_argument("-b", "--blockage_file", dest="blockage_file",
                        help="Masks for each transect identifying critical "
                             "land blockage.", metavar="BLOCKFILE",
                        required=True)
    args = parser.parse_args()

    dsMask = xarray.open_dataset(args.input_mask_filename)

    dsBlockages = xarray.open_dataset(args.blockage_file)

    dsMask = add_critical_land_blockages(dsMask, dsBlockages)
    dsMask.to_netcdf(args.output_mask_filename)
