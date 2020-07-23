from __future__ import print_function
import os
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma

#------------------------------------------------------------------------------------

def get_weights_filename(inputStrings):

    import hashlib

    hashInput = ""
    for inputString in inputStrings:
        hashInput = hashInput + inputString

    hashOutput = hashlib.md5(hashInput).hexdigest()
    weightFilename = "weights_%s.nc" %(hashOutput)

    return weightFilename

#------------------------------------------------------------------------------------

def generate_weights_file(SCRIPFilenameSrc, SCRIPFilenameDst, weightsFilename, reuseWeights):

    if (not reuseWeights or not os.path.isfile(weightsFilename)):

        cmd = "ESMF_RegridWeightGen --source %s --destination %s --weight %s --src_regional --dst_regional --ignore_unmapped" \
          %(SCRIPFilenameSrc, SCRIPFilenameDst, weightsFilename)
        print(cmd)
        os.system(cmd)

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

def regrid_mpas_array(mpasArray, weights):

    n_s      = weights["n_s"]

    col      = weights["col"]
    row      = weights["row"]
    S        = weights["S"]

    nColumns = weights["nColumns"]
    nRows    = weights["nRows"]

    # linear regrid
    mpasArrayRegridLinear = np.zeros(nRows*nColumns)

    for i_s in range(0,n_s):
        mpasArrayRegridLinear[row[i_s]-1] = mpasArrayRegridLinear[row[i_s]-1] + S[i_s] * mpasArray[col[i_s]-1]

    # get two dimensional grid
    mpasArrayRegrid = ma.zeros((nRows,nColumns))

    for iRow in range(0,nRows):
        for iColumn in range(0,nColumns):

            i = iColumn + iRow*nColumns

            mpasArrayRegrid[iRow,iColumn] = mpasArrayRegridLinear[i]

    return mpasArrayRegrid

#------------------------------------------------------------------------------------
