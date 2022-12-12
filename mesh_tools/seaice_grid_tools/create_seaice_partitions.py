#!/usr/bin/env python
import argparse
import sys

# parsing input
parser = argparse.ArgumentParser(description='Create sea-ice partitions.')

parser.add_argument('-m', '--outputmesh', dest="meshFilename",       required=True,  help='MPAS mesh file for destination regridding mesh.')
parser.add_argument('-o', '--outputDir',  dest="outputDir",          required=True,  help='Output directory for temporary files and partition files.')
parser.add_argument('-c', '--cullerDir',  dest="mpasCullerLocation", required=False, help='Location of MPAS MpasCellCuller.x executable.')
parser.add_argument('-p', '--prefix',     dest="outputPrefix",       required=False, help='prefix for output partition filenames.', default="graph.info")
parser.add_argument('-x', '--plotting',   dest="plotting",           required=False, help='create diagnostic plotting file of partitions', action='store_false')
parser.add_argument('-g', '--metis',      dest="metis",              required=False, help='name of metis utility', default="gpmetis")
parser.add_argument('-n', '--nProcs',     dest="nProcs",             required=False, help='number of processors to create partition for.', type=int)
parser.add_argument('-f', '--nProcsFile', dest="nProcsFile",         required=False, help='number of processors to create partition for.')

args = parser.parse_args()

# number of processors
nProcsArray = []
if (args.nProcs is not None and args.nProcsFile is None):
    nProcsArray.append(args.nProcs)
elif (args.nProcs is None and args.nProcsFile is not None):
    fileNProcs = open(args.nProcsFile,"r")
    nProcsLines = fileNProcs.readlines()
    fileNProcs.close()
    for line in nProcsLines:
        nProcsArray.append(int(line))
elif (args.nProcs is None and args.nProcsFile is None):
    print("ERROR: Must specify nProcs or nProcsFile")
    sys.exit()
elif (args.nProcs is not None and args.nProcsFile is not None):
    print("ERROR: Can't specify both nProcs or nProcsFile")
    sys.exit()

# create partitions
regionFilename = args.outputDir + "/regions.nc"
outputPrefix = args.outputDir + "/" + args.outputPrefix

from gen_seaice_mesh_partition import gen_seaice_mesh_partition
for nProcs in nProcsArray:
    gen_seaice_mesh_partition(args.meshFilename, regionFilename, nProcs, args.mpasCullerLocation, outputPrefix, args.plotting, args.metis, False)
