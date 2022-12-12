#!/usr/bin/env python
import argparse
import os, sys
import subprocess

# parsing input
parser = argparse.ArgumentParser(description='Perform preparatory work for making seaice partitions.')

parser.add_argument('-i', '--inputmesh',  dest="meshFilenameSrc", required=True, help='MPAS mesh file for source regridding mesh.')
parser.add_argument('-p', '--presence',   dest="filenameData",    required=True, help='Input ice presence file for source mesh.')
parser.add_argument('-m', '--outputmesh', dest="meshFilenameDst", required=True, help='MPAS mesh file for destination regridding mesh.')
parser.add_argument('-o', '--outputDir',  dest="outputDir",       required=True, help='Output directory for temporary files and partition files.')

args = parser.parse_args()

# Check if output directory exists
if (not os.path.isdir(args.outputDir)):
    print("ERROR: Output directory does not exist.")
    sys.exit()

# 1) Regrid the ice presence from the input data mesh to the grid of choice
from regrid_to_other_mesh import regrid_to_other_mesh
print("Regrid to desired mesh...")
filenameOut = args.outputDir + "/icePresent_regrid.nc"
regrid_to_other_mesh(args.meshFilenameSrc, args.filenameData, args.meshFilenameDst, filenameOut)


# 2) create icePresence variable
print("fix_regrid_output...")

# check executable exists
if (not os.path.exists("fix_regrid_output.exe")):
    print("ERROR: fix_regrid_output.exe does not exist.")
    sys.exit()

inputFile  = args.outputDir + "/icePresent_regrid.nc"
outputFile = args.outputDir + "/icePresent_regrid_modify.nc"
subprocess.call(["./fix_regrid_output.exe", inputFile, args.meshFilenameDst, outputFile])


# 3) create variable icePresenceExtended
from extend_seaice_mask import extend_seaice_mask
print("extend_seaice_mask...")
filenamePresence = args.outputDir + "/icePresent_regrid_modify.nc"
extend_seaice_mask(args.meshFilenameDst,filenamePresence,0.0,False)


# 4) Make the regions file from the icePresenceExtended variable
from make_regions_file import make_regions_file
print("make_regions_file...")
filenameIcePresent = args.outputDir + "/icePresent_regrid_modify.nc"
filenameOut = args.outputDir + "/regions.nc"
make_regions_file(filenameIcePresent, args.meshFilenameDst, "two_region_eq", "icePresenceExtended", 0.5, filenameOut)
