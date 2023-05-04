#!/usr/bin/env python

"""
This script copies a restart file of a MALI simulation
and re-calculates missing state variables for a 
missing time level and writes them to an updated restart file.
"""

import argparse
import os
import shutil
import xarray as xr
import numpy as np


def main():
    parser = argparse.ArgumentParser(
                        description='process MALI outputs for the ISMIP6'
                                    'submission')
    parser.add_argument("-f", "--file", dest="file_in",
                        required=True,
                        help="restart file to be read in")
    parser.add_argument("-o", "--output_file", dest="file_out",
                        required=True,
                        help="output file name")
    parser.add_argument("-p", "--output_file_path",
                        dest="output_path")
    
    args = parser.parse_args()
    
    # read in a restart file that needs to be re-written
    if args.file_in is None:
        print("--- restart file is not provided. Aborting... ---")
    else:
        print("\n--- Reading in the restart and output state files ---")

        file_in = xr.open_dataset(args.file_in, decode_times=False, decode_cf=False)

        # get needed info from restart file
        cellMask = file_in['cellMask'][:, :]
        thickness = file_in['thickness'][:,:]
        bedTopography = file_in['bedTopography'][:,:]
        sfcAirTemp = file_in['surfaceAirTemperature'][:,:]
        uReconstructX = file_in['uReconstructX'][:,:,:]
        uReconstructY = file_in['uReconstructY'][:,:,:]
        layerThicknessFractions = file_in['layerThicknessFractions']
        nTime = file_in.dims['Time']
        nCells = file_in.dims['nCells']
        nVertLevels = file_in.dims['nVertLevels']

        # xtime needs some massaging for xarray not to mangle it
        xtime = file_in['xtime']
        xtimeStr = xtime.data.tobytes().decode()  # convert to str
        xtime2 = xr.DataArray(np.array([xtimeStr], dtype = np.dtype(('S', 64))), dims = ['Time'])  # convert back to char array
        # followed example here: https://github.com/pydata/xarray/issues/3407

        floating_iceMask = (cellMask[:, :] & 4) // 4
        seaLevel = 0.0
        rhoi = 910.0
        rhoo = 1028.0

        print(f'nTime={nTime}, nCells={nCells}')

        layerInterfaceFractions = np.zeros(nVertLevels+1, dtype=float)
        lowerSfc = np.zeros([nTime, nCells], dtype=float)
        upperSfc = np.zeros([nTime, nCells], dtype=float)
        sfcTemp = np.zeros([nTime, nCells], dtype=float)
        xvelmean = np.zeros([nTime, nCells], dtype=float)
        yvelmean = np.zeros([nTime, nCells], dtype=float)
        # the following need to be in the file so ncrcat will work but processing won't use
        # values, so can leave as zeros
        surfaceSpeed = np.zeros([nTime, nCells], dtype=float)
        vonMisesStress = np.zeros([nTime, nCells], dtype=float)
        deltat = np.zeros([nTime,], dtype=float)
        daysSinceStart = np.zeros([nTime,], dtype=float)
        
        print("\n--- calculating the missing state variables ---")

        # layerInterfaceFractions are the fraction associated with each interface
        layerInterfaceFractions[0] = 0.5 * layerThicknessFractions[0]
        for k in range(1, nVertLevels):
            layerInterfaceFractions[k] = 0.5 * (layerThicknessFractions[k-1]
                                         + layerThicknessFractions[k])
        layerInterfaceFractions[nVertLevels] = 0.5 * layerThicknessFractions[nVertLevels-1]
        print("layerThicknessFractions:", layerThicknessFractions[:].data)
        print("layerInterfaceFractions:", layerInterfaceFractions)

        for i in range(nTime):
            # calculate surface temperature (unit in Kelvin)
            sfcTemp[i,:] = np.minimum(273.15, sfcAirTemp[i,:]) # 0 celsius = 273 Kelvin
            print('surfaceTemperature processed')

            lowerSfc[i,:] = np.where(floating_iceMask, seaLevel - thickness[i,:] * (rhoi / rhoo), bedTopography[i,:])
            upperSfc[i,:] = lowerSfc[i,:] + thickness[i,:]
            print('lower/upperSurface processed')

            xvelmean[i,:] = np.sum(uReconstructX[i,:,:] * layerInterfaceFractions[:], axis=1)
            yvelmean[i,:] = np.sum(uReconstructY[i,:,:] * layerInterfaceFractions[:], axis=1)
            print('x/yvelmean processed')

        # create variable dictionary of fields to include in the new file
        # Note: ncrcat does not require that time-independent fields be in both
        # files, so we don't need to include them in the new file.
        out_data_vars = {
            'lowerSurface': (['Time', 'nCells'], lowerSfc),
            'upperSurface': (['Time', 'nCells'], upperSfc),
            'surfaceTemperature': (['Time', 'nCells'], sfcTemp),
            'xvelmean': (['Time', 'nCells'], xvelmean),
            'yvelmean': (['Time', 'nCells'], yvelmean),
            'surfaceSpeed': (['Time', 'nCells'], surfaceSpeed),
            'vonMisesStress': (['Time', 'nCells'], vonMisesStress),
            'deltat': (['Time',], deltat ),
            'daysSinceStart': (['Time',], daysSinceStart),
            'xtime': xtime2
            }
        dataOut = xr.Dataset(data_vars=out_data_vars)  # create xarray dataset object

        print("\n--- copying over dimensions from the restart file ---")
        for dim in file_in.dims:
            if not dim in dataOut.dims:
                dataOut.expand_dims({dim: file_in.dims[dim]})
                print('new file dimension', dataOut.dims)
                print("   Copying dimension", dim)

        dataOut.xtime.encoding.update({"char_dim_name": "StrLen"})  # another hacky thing to make xarray handle xtime correctly
        # learned this from: https://github.com/pydata/xarray/issues/2895

        print("\n--- copying over unmodified variables from the restart file ---")
        for var in ['thickness', 'uReconstructX', 'uReconstructY', 'bedTopography',
                    'basalTemperature', 'betaSolve', 'cellMask', 'damage']:
            print("   Copying", var)
            dataOut[var] = file_in[var]

        # save/write out the new file
        # define the path to which the output (processed) files will be saved
        if args.output_path is None:
            output_path = os.getcwd()
        else:
            output_path = args.output_path
        
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        print(f"file output path: {output_path}")
        file_out_path = os.path.join(output_path, args.file_out)
        dataOut.to_netcdf(file_out_path, mode='w', unlimited_dims=['Time'])
        file_in.close()
        
        print("\n--- process complete! ---")

if __name__ == "__main__":
    main()
