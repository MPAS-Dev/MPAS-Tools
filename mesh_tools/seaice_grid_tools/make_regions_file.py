#!/usr/bin/env python
import argparse

from seaice.regions import make_regions_file


if __name__ == "__main__":

    # input parsing
    parser = argparse.ArgumentParser(description='Make a partition regions file.')

    parser.add_argument('-p', '--presence', dest="filenameIcePresent", required=True, help='ice presence file')
    parser.add_argument('-m', '--mesh',     dest="filenameMesh",       required=True, help='MPAS mesh file')
    # region type options: "two_region_eq", "three_region", "three_region_eq", "five_region_eq"
    parser.add_argument('-t', '--type',     dest="regionType",         required=True, help='region type')
    parser.add_argument('-v', '--varname',  dest="varname",            required=True, help='presence var name')
    parser.add_argument('-o', '--output',   dest="filenameOut",        required=True, help='output regions file')
    parser.add_argument('-l', '--limit',    dest="limit",              required=False, default=0.5, type=float, help='presence region limit')

    args = parser.parse_args()

    make_regions_file(args.filenameIcePresent, args.filenameMesh, args.regionType, args.varname, args.limit, args.filenameOut)
