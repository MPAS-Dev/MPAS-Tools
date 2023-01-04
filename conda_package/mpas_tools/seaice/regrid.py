import os
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import subprocess
import hashlib

from .mesh import make_mpas_scripfile_on_cells


#------------------------------------------------------------------------------------

def get_weights_filename(inputStrings):

    hashInput = ""
    for inputString in inputStrings:
        hashInput = hashInput + inputString

    hashOutput = hashlib.md5(hashInput).hexdigest()
    weightFilename = "weights_%s.nc" %(hashOutput)

    return weightFilename

#------------------------------------------------------------------------------------

def generate_weights_file(SCRIPFilenameSrc, SCRIPFilenameDst, weightsFilename, reuseWeights):

    if not reuseWeights or not os.path.isfile(weightsFilename):

        args = ["ESMF_RegridWeightGen",
                "--source", SCRIPFilenameSrc,
                "--destination", SCRIPFilenameDst,
                "--weight", weightsFilename,
                "--src_regional", "--dst_regional", "--ignore_unmapped"]

        subprocess.run(args, check=True)

#------------------------------------------------------------------------------------

def load_weights_file(weightsFilename):

    weights = {}

    weightsFile = Dataset(weightsFilename,"r")

    weights["n_s"] = len(weightsFile.dimensions["n_s"])

    weights["col"] = weightsFile.variables["col"][:]
    weights["row"] = weightsFile.variables["row"][:]
    weights["S"]   = weightsFile.variables["S"][:]

    dst_grid_dims = weightsFile.variables["dst_grid_dims"][:]
    weights["nRows"]    = dst_grid_dims[0]
    weights["nColumns"] = dst_grid_dims[1]

    weightsFile.close()

    return weights

#------------------------------------------------------------------------------------

def regrid_obs_array(obsArray, weights):

    n_s      = weights["n_s"]

    col      = weights["col"]
    row      = weights["row"]
    S        = weights["S"]

    nColumns = weights["nColumns"]
    nRows    = weights["nRows"]

    nRowsIn    = obsArray.shape[0]
    nColumnsIn = obsArray.shape[1]

    # get two dimensional grid
    obsArrayRegrid = ma.zeros((nRows,nColumns))

    for i_s in range(0,n_s):

        iRow      = (row[i_s]-1) / nRows
        iColumn   = (row[i_s]-1) % nRows

        iRowIn    = (col[i_s]-1) / nRowsIn
        iColumnIn = (col[i_s]-1) % nRowsIn

        obsArrayRegrid[iRow,iColumn] = obsArrayRegrid[iRow,iColumn] + S[i_s] * obsArray[iRowIn,iColumnIn]

    return obsArrayRegrid

#------------------------------------------------------------------------------------

def regrid_mpas_array(weightsFilename, mpasArrayIn):

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

#---------------------------------------------------------------------------------

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
    iceFractionOut, iceFractionOutMask = regrid_mpas_array(weightsFilename, iceFractionIn)



    # output
    print("Output...")
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells",nCells)

    iceFractionVar = fileOut.createVariable("iceFraction","d",dimensions=["nCells"])
    iceFractionVar[:] = iceFractionOut[:]

    iceFractionMaskVar = fileOut.createVariable("iceFractionMask","i",dimensions=["nCells"])
    iceFractionMaskVar[:] = iceFractionOutMask[:]

    fileOut.close()
