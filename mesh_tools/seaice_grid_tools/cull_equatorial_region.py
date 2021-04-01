from __future__ import print_function
import argparse
import os, sys
import subprocess
import shutil
from netCDF4 import Dataset
import numpy as np

# parsing input
parser = argparse.ArgumentParser(description='Perform prepatory work for making seaice partitions.')

parser.add_argument('-m', '--presenceMesh', dest="presenceMesh",       required=True,  help='MPAS mesh file for source regridding mesh.')
parser.add_argument('-p', '--presence',     dest="presenceDatafile",   required=True,  help='Input ice presence file for source mesh.')
parser.add_argument('-t', '--tmpDir',       dest="tmpDir",             required=True,  help='Tmp directory for temporary files.')
parser.add_argument('-c', '--culler',       dest="mpasCullerLocation", required=False, help='location of cell culler')
parser.add_argument('-e', '--extendDist',   dest="extendDistance",     required=True,  help='distance (km) to extend ice present region', type=float)
parser.add_argument('-i', '--inputMesh',    dest="meshFilenameSrc",    required=True,  help='Input mesh filename to cull.')
parser.add_argument('-o', '--outputMesh',   dest="meshFilenameDst",    required=True,  help='Culled output mesh filename.')

args = parser.parse_args()

# Check if output directory exists
if (not os.path.isdir(args.tmpDir)):
    print("ERROR: Tmp directory does not exist: ", args.tmpDir)
    sys.exit()

# arguments
if (args.mpasCullerLocation == None):
    meshToolsDir = os.path.dirname(os.path.realpath(__file__)) + "/../mesh_conversion_tools/"
else:
    meshToolsDir = args.mpasCullerLocation
if (not os.path.exists(meshToolsDir + "/MpasCellCuller.x")):
    print("ERROR: MpasCellCuller.x does not exist at the requested loaction.")
    sys.exit();


# 1) Regrid the ice presence from the input data mesh to the grid of choice
from regrid_to_other_mesh import regrid_to_other_mesh
print("Regrid to desired mesh...")
filenameOut = args.tmpDir + "/icePresent_regrid.nc"
regrid_to_other_mesh(args.presenceMesh, args.presenceDatafile, args.meshFilenameSrc, filenameOut)


# 2) create icePresence variable
print("fix_regrid_output...")

# check executable exists
if (not os.path.exists("fix_regrid_output.exe")):
    print("ERROR: fix_regrid_output.exe does not exist.")
    sys.exit()

inputFile  = args.tmpDir + "/icePresent_regrid.nc"
outputFile = args.tmpDir + "/icePresent_regrid_modify.nc"
subprocess.call(["./fix_regrid_output.exe", inputFile, args.meshFilenameSrc, outputFile])


# 3) create variable icePresenceExtended
from extend_seaice_mask import extend_seaice_mask
print("extend_seaice_mask...")
filenamePresence = args.tmpDir + "/icePresent_regrid_modify.nc"
extend_seaice_mask(args.meshFilenameSrc,filenamePresence,args.extendDistance,False)


# 4) copy input mesh file to output file
shutil.copyfile(args.meshFilenameSrc, args.meshFilenameDst)


# 5) add cullCell variable
fileIn = Dataset(filenamePresence,"r")
icePresenceExtended = fileIn.variables["icePresenceExtended"][:]
fileIn.close()

fileIn = Dataset(args.meshFilenameSrc,"a")

nCells = len(fileIn.dimensions["nCells"])
cullCell = np.zeros(nCells)
for iCell in range(0,nCells):
    cullCell[iCell] = 1 - icePresenceExtended[iCell]

try:
    cullCellVar = fileIn.createVariable("cullCell","i4",dimensions=["nCells"])
except:
    cullCellVar = fileIn.variables["cullCell"]
cullCellVar[:] = cullCell[:]

fileIn.close()


# 6) cull equatorial cells
os.system("%s/MpasCellCuller.x %s %s -c" %(meshToolsDir,args.meshFilenameSrc,args.meshFilenameDst))
