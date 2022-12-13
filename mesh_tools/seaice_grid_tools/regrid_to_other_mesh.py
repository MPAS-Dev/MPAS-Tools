#!/usr/bin/env python
import argparse

from mpas_tools.seaice.regrid import regrid_to_other_mesh


if __name__ == "__main__":

    # parsing input
    parser = argparse.ArgumentParser(description='Regrid sea ice presence to other meshes.')

    parser.add_argument('-i', '--inputmesh',  dest="meshFilenameSrc", required=True, help='MPAS mesh file for source regridding mesh')
    parser.add_argument('-p', '--presence',   dest="filenameData",    required=True, help='Input ice presence file for source mesh')
    parser.add_argument('-m', '--outputmesh', dest="meshFilenameDst", required=True, help='MPAS mesh file for destination regridding mesh')
    parser.add_argument('-o', '--output',     dest="filenameOut",     required=True, help='Output ice presence file for destination mesh')

    args = parser.parse_args()

    regrid_to_other_mesh(args.meshFilenameSrc, args.filenameData, args.meshFilenameDst, args.filenameOut)
