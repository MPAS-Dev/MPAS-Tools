#!/usr/bin/env python
import argparse

from seaice.partition import gen_seaice_mesh_partition


if __name__ == "__main__":

    # parsing
    parser = argparse.ArgumentParser(description='Create sea ice grid partition')

    parser.add_argument('-m', '--mesh',        dest="meshFilename",         required=True,  help='MPAS mesh file')
    parser.add_argument('-r', '--regions',     dest="regionFilename",       required=True,  help='region file')
    parser.add_argument('-n', '--nprocs',      dest="nProcs",               required=True,  help='number of processors', type=int)
    parser.add_argument('-c', '--culler',      dest="mpasCullerLocation",   required=False, help='location of cell culler')
    parser.add_argument('-o', '--outprefix',   dest="outputPrefix",         required=False, help='output graph file prefic', default="graph.info")
    parser.add_argument('-p', '--plotting',    dest="plotting",             required=False, help='create diagnostic plotting file of partitions', action='store_true')
    parser.add_argument('-g', '--metis',       dest="metis",                required=False, help='name of metis utility', default="gpmetis")
    parser.add_argument('-e', '--equatorcull', dest="cullEquatorialRegion", required=False, help='create diagnostic plotting file of partitions', action='store_true')

    args = parser.parse_args()

    gen_seaice_mesh_partition(args.meshFilename, args.regionFilename, args.nProcs, args.mpasCullerLocation, args.outputPrefix, args.plotting, args.metis, args.cullEquatorialRegion)
