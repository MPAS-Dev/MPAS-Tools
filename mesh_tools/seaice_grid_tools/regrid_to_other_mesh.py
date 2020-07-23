from __future__ import print_function
from netCDF4 import Dataset
import numpy as np
import sys
import argparse
import os

from seaice_regrid import *
from seaice_mesh import *

#--------------------------------------------------------------

def regrid_mpas_array(weightsFilename, mpasArrayIn, nCells, nEdgesOnCell, cellsOnCell, fixPolarCoasts):

    nMax = 100

    # load weights
    print("Load weights...", weightsFilename)
    weightsFile = Dataset(weightsFilename,"r")

    n_s = len(weightsFile.dimensions["n_s"])
    n_a = len(weightsFile.dimensions["n_a"])
    n_b = len(weightsFile.dimensions["n_b"])

    col = weightsFile.variables["col"][:]
    row = weightsFile.variables["row"][:]
    S   = weightsFile.variables["S"][:]

    weightsFile.close()

    # regrid
    print("Regrid array...")
    mpasArrayOut = np.zeros(n_b)
    mpasArrayOutMask = np.zeros(n_b,dtype="i")
    for i_s in range(0,n_s):
        mpasArrayOut[row[i_s]-1] = mpasArrayOut[row[i_s]-1] + S[i_s] * mpasArrayIn[col[i_s]-1]
        mpasArrayOutMask[row[i_s]-1] = 1

    return mpasArrayOut, mpasArrayOutMask

#--------------------------------------------------------------

def regrid_to_other_mesh(meshFilenameSrc, filenameData, meshFilenameDst, filenameOut):

    # make scrip files
    print("Make scrip files...")
    SCRIPFilenameSrc = "scrip_src_tmp.nc"
    SCRIPFilenameDst = "scrip_dst_tmp.nc"

    titleSrc = "MPAS grid src"
    titleDst = "MPAS grid dst"

    make_mpas_scripfile_on_cells(meshFilenameSrc, SCRIPFilenameSrc, titleSrc)
    make_mpas_scripfile_on_cells(meshFilenameDst, SCRIPFilenameDst, titleDst)

    # generate weights file
    print("Generate weights...")
    weightsFilename = os.getcwd() + "/weights_tmp.nc"
    generate_weights_file(SCRIPFilenameSrc, SCRIPFilenameDst, weightsFilename, False)


    # load output mesh
    print("Load output mesh...")
    meshFile = Dataset(meshFilenameDst,"r")

    nCells = len(meshFile.dimensions["nCells"])

    nEdgesOnCell = meshFile.variables["nEdgesOnCell"][:]
    cellsOnCell = meshFile.variables["cellsOnCell"][:]
    cellsOnCell[:] = cellsOnCell[:] - 1

    latCell = meshFile.variables["latCell"][:]
    lonCell = meshFile.variables["lonCell"][:]

    meshFile.close()

    # load data
    print("Load input data...")
    fileIn = Dataset(filenameData,"r")

    iceFractionIn = fileIn.variables["iceFraction"][:]

    fileIn.close()



    # regrid
    print("Regrid array...")
    iceFractionOut, iceFractionOutMask = regrid_mpas_array(weightsFilename, iceFractionIn, nCells, nEdgesOnCell, cellsOnCell, False)



    # output
    print("Output...")
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells",nCells)

    iceFractionVar = fileOut.createVariable("iceFraction","d",dimensions=["nCells"])
    iceFractionVar[:] = iceFractionOut[:]

    iceFractionMaskVar = fileOut.createVariable("iceFractionMask","i",dimensions=["nCells"])
    iceFractionMaskVar[:] = iceFractionOutMask[:]

    fileOut.close()

#--------------------------------------------------------------

if __name__ == "__main__":

    # parsing input
    parser = argparse.ArgumentParser(description='Regrid sea ice presence to other meshes.')

    parser.add_argument('-i', '--inputmesh',  dest="meshFilenameSrc", required=True, help='MPAS mesh file for source regridding mesh')
    parser.add_argument('-p', '--presence',   dest="filenameData",    required=True, help='Input ice presence file for source mesh')
    parser.add_argument('-m', '--outputmesh', dest="meshFilenameDst", required=True, help='MPAS mesh file for destination regridding mesh')
    parser.add_argument('-o', '--output',     dest="filenameOut",     required=True, help='Output ice presence file for destination mesh')

    args = parser.parse_args()

    regrid_to_other_mesh(args.meshFilenameSrc, args.filenameData, args.meshFilenameDst, args.filenameOut)
