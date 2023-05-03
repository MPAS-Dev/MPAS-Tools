#!/usr/bin/env python

"""
This script copies a restart file of a MALI simulation
and re-calculates missing state variables for a 
missing time level and writes them to an updated restart file.
"""

import argparse
from subprocess import check_call
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
        print("\n--- Reading in/copying the restart file ---")

        file_in = xr.open_dataset(args.file_in, engine="netcdf4", decode_cf=False)
        file_in_copy = file_in.copy(deep=True)
        
        # calculate mising fields for the missing time level
        cellMask = file_in['cellMask'][:, :]
        thickness = file_in['thickness'][:,:]
        bedTopography = file_in['thickness'][:,:]
        sfcAirTemp = file_in['surfaceAirTemperature'][:,:]
        uReconstructX = file_in['uReconstructX'][:,:]
        uReconstructY = file_in['uReconstructY'][:,:]
        layerThicknessFractions = file_in['layerThicknessFractions']
        nTime = file_in.dims['Time']
        nCells = file_in.dims['nCells']
        nVertLevels = file_in.dims['nVertLevels']
        floating_iceMask = (cellMask[:, :] & 4) / 4 * (cellMask[:, :] & 2) / 2  # floating and dynamic
        seaLevel = 0.0
        rhoi = 910.0
        rhoo = 1028.0
        
        layerInterfaceFractions = np.zeros(nVertLevels+1, dtype=float)
        lowerSfc = np.zeros([nTime, nCells], dtype=float)
        upperSfc = np.zeros([nTime, nCells], dtype=float) 
        sfcTemp = np.zeros([nTime, nCells], dtype=float) 
        xvelmean = np.zeros([nTime, nCells], dtype=float) 
        yvelmean = np.zeros([nTime, nCells], dtype=float) 
        
        print("\n--- calculating the missing state variables ---")

        # layerInterfaceFractions are the fraction associated with each interface
        layerInterfaceFractions[0] = 0.5 * layerThicknessFractions[0]
        for k in range(1, nVertLevels):
            layerInterfaceFractions[k] = 0.5 * (layerThicknessFractions[k-1]
                                         + layerThicknessFractions[k])
        layerInterfaceFractions[nVertLevels] = 0.5 * layerThicknessFractions[nVertLevels-1]
        
        for i in range(nTime):
            for j in range(nCells):
                # calculate lower and upper surface elevation
                if (cellMask[i,j] & 4) // 4 == 1: # floating ice
                    lowerSfc[i,j] = seaLevel - thickness[i,j] * (rhoi / rhoo)
                else: 
                    lowerSfc[i,j] = bedTopography[i,j]
            
                upperSfc[i,j] = lowerSfc[i,j] + thickness[i,j]
             
                # calculate surface temperature (unit in Kelvin)
                sfcTemp[i,j] = np.minimum(273.15, sfcAirTemp[i,j]) # 0 celsius = 273 Kelvin

                # calculate vertical mean velocity
                xvelmean[i,j] = np.sum(uReconstructX[i,j] * layerInterfaceFractions[:])
                yvelmean[i,j] = np.sum(uReconstructY[i,j] * layerInterfaceFractions[:])
        
        file_in_copy['lowerSurface'] = (['Time', 'nCells'], lowerSfc)
        file_in_copy['upperSurface'] = (['Time', 'nCells'], upperSfc)
        file_in_copy['surfaceTemperature'] = (['Time', 'nCells'], sfcTemp)
        file_in_copy['xvelmean'] = (['Time', 'nCells'], xvelmean)
        file_in_copy['yvelmean'] = (['Time', 'nCells'], yvelmean)
        
        print("\n--- writing out to a new file ---")
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
        file_in_copy.to_netcdf(file_out_path, mode='a')
        file_in.close()
        
        print("\n--- process compleate! ---")

if __name__ == "__main__":
    main()