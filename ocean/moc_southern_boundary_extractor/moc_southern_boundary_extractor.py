#!/usr/bin/env python

'''
This script takes a mesh file (-m flag) and a file with MOC regions masks
(-f flag) produce by the MPAS mask creator.  The script produces a copy of
the contents of the MOC mask file, adding transects that mark the southern
boundary of each region in a file indicated with the -o flag.  The transect
is applied only to vertices and edges, not cells, because the need for southern
boundary transect data on cells is not foreseen.

Author: Xylar Asay-Davis
last modified: 5/22/2018
'''

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import argparse

from mpas_tools.io import \
    write_netcdf
from mpas_tools.ocean.moc import \
    add_moc_southern_boundary_transects


if __name__ == "__main__":

    parser = \
        argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f', '--in_file', dest='in_file',
                        help='Input file with MOC masks', metavar='IN_FILE',
                        required=True)
    parser.add_argument('-m', '--mesh_file', dest='mesh_file',
                        help='Input mesh file', metavar='MESH_FILE',
                        required=True)
    parser.add_argument('-o', '--out_file', dest='out_file',
                        help='Output file for MOC masks and southern-boundary '
                        'transects', metavar='OUT_FILE',
                        required=True)
    args = parser.parse_args()

    dsMasks = xarray.open_dataset(args.in_file)
    dsMesh = xarray.open_dataset(args.mesh_file)

    dsMasksAndTransects = add_moc_southern_boundary_transects(dsMasks, dsMesh)

    write_netcdf(dsMasksAndTransects, args.out_file)
