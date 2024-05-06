#!/usr/bin/env python

import os
import shutil
import xarray as xr
import subprocess

from mpas_tools.logging import check_call
from mpas_tools.scrip.from_mpas import scrip_from_mpas
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--oceanMesh", dest="mpasMeshFile", help="MPAS-Ocean mesh to be interpolated from", metavar="FILENAME")
parser.add_option("--oceanDiag", dest="mpasDiagFile", help="MPAS-Ocean diagnostic output file", metavar="FILENAME")
parser.add_option("--ice", dest="maliFile", help="MALI mesh to be interpolated onto", metavar="FILENAME")
parser.add_option("-n", "--ntasks", dest="ntasks", help="Number of processors to use with ESMF_RegridWeightGen")
parser.add_option("-m", "--method", dest="method", help="Remapping method, either 'bilinear' or 'neareststod'")
parser.add_option("-o", "--outFile",dest="outputFile", help="Desired name of output file", metavar="FILENAME", default="mpas_to_mali_remapped.nc")
options, args = parser.parse_args()

#make copy of mpasMeshFile
tmp_mpasMeshFile = "tmp_mpasMeshFile.nc"
subprocess.run(["ncks", options.mpasMeshFile, tmp_mpasMeshFile]) 

#Preprocessing:
#TO DO: Create capability to read multiple files and time average over relavant fields

#TO DO: Vertically interpolate between MPAS-O to MALI grids

#main:
#calculate thermal forcing
calc_ocean_thermal_forcing(options.mpasMeshFile, options.mpasDiagFile)

#TO DO: copy landIceFreshwaterFlux to mpas mesh file copy for interpolation. Make name relevant to MALI

#remap mpas to mali
remap_mpas_to_mali(options.mpasMeshFile, options.maliMeshFile, options.ntasks, options.method, options.outputFile)

def calc_ocean_thermal_forcing(tmp_mpasMeshFile, mpasDiagFile)
    print("Calculating Thermal Forcing")
    print("... gathering data ...")
    fM = xr.open_dataset(options.mpasMeshFile, decode_times=False, decode_cf=False)
    minLevelCell = fM['minLevelCell'][:].data
    maxLevelCell = fM['maxLevelCell'][:].data
    layerThickness = fM['layerThickness'][:,:,:].data 
    landIceFloatingMask = fM['landIceFloatingMask'][:,:].data
    nCells = fM.sizes['nCells']

    fD = xr.open_dataset(options.mpasDiagFile, decode_times=False, decode_cf=False)
    temperature = fD['timeMonthly_avg_activeTracers_temperature'][:,:,:].data
    salinity = fD['timeMonthly_avg_activeTracers_salinity'][:,:,:].data
    density = fD['time_monthy_avg_density'][:,:,:].data
    atmPressure = fD['time_monthly_avg_atmosphericPressure'][:,:].data
    coeff_0_openOcean = fD.attrs['config_open_ocean_freezing_temperature_coeff_0']
    coeff_S_openOcean = fD.attrs['config_open_ocean_freezing_temperature_coeff_S']
    coeff_p_openOcean = fD.attrs['config_open_ocean_freezing_temperature_coeff_p']
    coeff_pS_openOcean = fD.attrs['config_open_ocean_freezing_temperature_coeff_pS']
    coeff_mushy_az1_liq = fD.attrs['config_open_ocean_freezing_temperature_coeff_mushy_az1_liq']
    coeff_mushy_openOcean = 1/coeff_mushy_az1_liq 
    coeff_0_cavity = fD.attrs['config_land_ice_cavity_freezing_temperature_coeff_0']
    coeff_S_cavity = fD.attrs['config_land_ice_cavity_freezing_temperature_coeff_S']
    coeff_p_cavity = fD.attrs['config_land_ice_cavity_temperature_coeff_p']
    coeff_pS_cavity = fD.attrs['config_land_ice_cavity_freezing_temperature_coeff_pS']
    coeff_mushy_az1_liq = fD.attrs['config_land_ice_cavity_freezing_temperature_coeff_mushy_az1_liq']
    coeff_mushy_cavity = 1/coeff_mushy_az1_liq
    Time = fD.sizes['Time']

    gravity = 9.81

    for iTime in range(Time-1):

        #calculate pressure: NEED TO CHECK KMIN AND K INDEXING
        for iCell in range(nCells-1):
                kmin = minLevelCell[iCell]
                kmax = maxLevelCell[iCell]
                pressure[iTime,kmin,iCell] = atmPressure[iTime,iCell] + density[iTime,iCell,kmin]*gravity*0.5*layerThickness[iTime,iCell,kmin]

                for k in np.arange(kmin, kmax):
                    pressure[iTime,k,iCell] = pressure[iTime,iCell,k] + 0.5*gravity*(density[iTime,iCell,k]*layerThickness[iTime,iCell,k] + density[iTime,iCell,k-1]*layerThickness[iTime,iCell,k-1])

                if (landIceFloatingMask[iTime,iCell] == 1):
                    ocn_freezing_temperature[iTime,iCell,:] = coeff_0_openOcean + coeff_S_openOcean * salinity[iTime,iCell,:] + coeff_p_openOcean * pressure[iTime,iCell,:] \
                            + coeff_pS_openOcean * pressure[iTime,iCell,:] * salinity[iTime,iCell,:] + coeff_mushy_openOcean * salinity[iTime,iCell,:] / (1.0 - salinity[iTime,iCell,:] / 1e3)
                
                elif (landIceFloatingMask[iTime,iCell] == 0):
                    ocn_freezing_temperature[iTime,iCell,:] = coeff_0_cavity + coeff_S_cavity * salinity[iTime,iCell,:] + coeff_p_cavity * pressure[iTime,iCell,:] \
                            + coeff_pS_cavity * pressure[iTime,iCell,:] * salinity[iTime,iCell,:] + coeff_mushy_cavity * salinity[iTime,iCell,:] / (1.0 - salinity[iTime,iCell,:] / 1e3)

    # Calculate thermal forcing
    oceanThermalForcing = temperature - ocn_freezing_temperature
    
    # Save thermal forcing to mesh file
    #TO DO: Change name to variable that MALI will recognize
    tf = xr.DataArray(oceanThermalForcing, dims=('Time','nCells','nVertLevels') 
    fM['oceanThermalForcing'] = tf
    fM.to_netcdf(tmp_mpasMeshFile)
    fM.close()

def remap_mpas_to_mali(mpas_meshFile, mali_meshFile, ntasks, method, outputFile)
    #create scrip files
    print("Creating Scrip Files")
    mali_scripfile = "tmp_mali_scrip.nc"
    mpas_scripfile = "tmp_mpas_scrip.nc"

    scrip_from_mpas(mali_meshFile, mali_scripfile)
    scrip_from_mpas(mpas_meshFile, mpas_scripfile)

    #create mapping file
    print("Creating Mapping File")
    args_esmf = ['srun', 
            '-n', ntasks, 'ESMF_RegridWeightGen',
            '--source', mpas_scripfile,
            '--destination', mali_scripfile,
            '--weight', "tmp_mapping_file.nc",
            '--method', method,
            '-i', '-64bit_offset',
            "--dst_regional", "--src_regional", '--ignore_unmapped']
    check_call(args_esmf)

    #TO DO: Create mask indentifying valid overlapping ocean cells

    #permute to work with ncremap
    shutil.copy(mpas_meshFile, "tmp_" + mpas_meshFile)
    tmp_mpasFile = "tmp_" + mpas_meshFile

    subprocess.run(["ncpdq", "-O", "-a", "Time,nVertLevels,nCells", tmp_mpasFile, tmp_mpasFile])
    subprocess.run(["ncpdq", "-O", "-a", "Time,nVertLevelsP1,nCells", tmp_mpasFile, tmp_mpasFile])
    subprocess.run(["ncpdq", "-O", "-a", "Time,nVertInterfaces,nCells", tmp_mpasFile, tmp_mpasFile])

    # remap the input data
    args_remap = ["ncremap",
            "-i", tmp_mpasFile,
            "-o", outputFile,
            "-m", "tmp_mapping_file.nc"]
    check_call(args_remap)

    #remove temporary files
    os.remove(mpas_scripfile)
    os.remove(mali_scripfile)
    os.remove(tmp_mpasFile)
    os.remove("tmp_mapping_file.nc")

