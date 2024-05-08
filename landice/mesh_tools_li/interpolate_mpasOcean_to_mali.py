#!/usr/bin/env python

import os
import shutil
import xarray as xr
import subprocess
import numpy as np
from scipy import stats as st
from mpas_tools.logging import check_call
from mpas_tools.scrip.from_mpas import scrip_from_mpas
from optparse import OptionParser

class mpasToMaliInterp:
    
    def __init__(self):
        parser = OptionParser()
        parser.add_option("--oceanMesh", dest="mpasMeshFile", help="MPAS-Ocean mesh to be interpolated from", metavar="FILENAME")
        parser.add_option("--oceanDiags", dest="mpasDiagsDir", help="Directory path where MPAS-Ocean diagnostic files are stored. Should be only netcdf files in that directory", metavar="FILENAME")
        parser.add_option("--ice", dest="maliFile", help="MALI mesh to be interpolated onto", metavar="FILENAME")
        parser.add_option("-n", "--ntasks", dest="ntasks", help="Number of processors to use with ESMF_RegridWeightGen")
        parser.add_option("-m", "--method", dest="method", help="Remapping method, either 'bilinear' or 'neareststod'")
        parser.add_option("-o", "--outFile",dest="outputFile", help="Desired name of output file", metavar="FILENAME", default="mpas_to_mali_remapped.nc")
        parser.add_option("-a","--yearlyAvg", dest="yearlyAvg", help="true or false option to average monthly output to yearly intervals", default="true")
        self.options, args = parser.parse_args()
        # << TO DO >>: Get rid of options sub-class, define options directly in self class

        # open and concatenate diagnostic dataset
        DS = xr.open_mfdataset(self.options.mpasDiagsDir + '/' + '*.nc', combine='nested', concat_dim='Time', decode_timedelta=False)        
        #open mpas ocean mesh
        OM = xr.open_dataset(tmp_mpasMeshFile, decode_times=False, decode_cf=False)
        #open MALI mesh
        IM = xr.open_dataset(options.maliMeshFile, decode_times=False, decode_cf=False)
        
        # variables for time averaging
        # << TO DO >>: Adjust for variables with no time dimension after open_mfdataset (e.g., Matt's areaCell example)
        self.temperature = DS['timeMonthly_avg_activeTracers_temperature'][:,:,:].data
        self.salinity = DS['timeMonthly_avg_activeTracers_temperature'][:,:,:].data
        self.density = DS['timeMonthly_avg_density'][:,:,:].data
        self.daysSinceStart = DS['daysSinceStart'][:].data
        self.stTime = DS.attrs['config_start_time']
        self.landIceFloatingMask = OM['landIceFloatingMiask'][:,:,:].data
        self.layerThickness = OM['layerThickness'][:,:].data 
        # << TO DO >> : add landIceFreshwaterFlux field
        # << TO DO >>: Figure out how to import xtime with xarray

        # additional variables for computing thermal forcing
        self.minLevelCell = OM['minLevelCell'][:].data
        self.maxLevelCell = OM['maxLevelCell'][:].data
        self.nCells = OM.sizes['nCells']
        self.coeff_0_openOcean = DS.attrs['config_open_ocean_freezing_temperature_coeff_0']
        self.coeff_S_openOcean = DS.attrs['config_open_ocean_freezing_temperature_coeff_S']
        self.coeff_p_openOcean = DS.attrs['config_open_ocean_freezing_temperature_coeff_p']
        self.coeff_pS_openOcean = DS.attrs['config_open_ocean_freezing_temperature_coeff_pS']
        self.coeff_mushy_az1_liq = DS.attrs['config_open_ocean_freezing_temperature_coeff_mushy_az1_liq']
        self.coeff_mushy_openOcean = 1/coeff_mushy_az1_liq 
        self.coeff_0_cavity = DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_0']
        self.coeff_S_cavity = DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_S']
        self.coeff_p_cavity = DS.attrs['config_land_ice_cavity_temperature_coeff_p']
        self.coeff_pS_cavity = DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_pS']
        self.coeff_mushy_az1_liq = DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_mushy_az1_liq']
        self.coeff_mushy_cavity = 1/coeff_mushy_az1_liq
        self.Time = DS.sizes['Time']

    def time_average_output(self)
        if (self.options.yearlyAvg == 'true'):
            yearsSinceStart = self.daysSinceStart / 365.0
            finalYear = np.floor(np.max(yearsSinceStart))
            timeStride = 1 # 1 year average

            years = np.arange(0, finalYear, timeStride)
            for i in range(len(years)-1):
                ind = np.where(yearsSinceStart >= years(i) and yearsSinceStart < years(i+1))
                self.newTemp[i,:,:] = np.mean(self.temperature[ind,:,:], axis=0)
                self.newSal[i,:,:] = np.mean(self.salinity[ind,:,:], axis=0)
                self.newDens[i,:,:] = np.mean(self.density[ind,:,:], axis=0)
                self.newLThick[i,:,:] = np.mean(self.layerThickness[ind,:,:], axis=0)
                self.newAtmPr[i,:] = np.mean(self.atmPressure[ind,:], axis=0)
                self.newFloatIceMask[i,:] = st.mode(self.landIceFloatingMask[ind,:], axis=0)

            # Build new xtime string. Defined on first day of each year
            stTime = np.datetime64(self.stTime[0:4])
            stYear = np.datetime_as_string(stTime,unit='Y').astype('float64')
            yearVec = yearsSinceStart + stYear
            self.newXtime = np.array([f"{int(year)}-01-01_00:00:00" for year in yearVec],dtype=np.str_)
        else #do nothing if not time averaging
            self.newTemp = self.temperature
            self.newSal = self.salinity
            self.newDens= self.density
            self.newLThick = self.layerThickness
            self.newAtmPr = self.atmPressure
            self.newFloatIceMask = self.landIceFloatingMask
            # << TO DO >>: Figure out what to do with xtime if no change. How do we get xarray to load original xtime variable?            
    
    def calc_ocean_thermal_forcing(self)

        gravity = 9.81

        for iTime in range(self.Time-1): # << TO DO >>: "Time" here needs to change to be consistent with time averaging

            #calculate pressure: 
            # << TO DO >>: NEED TO CHECK KMIN AND K INDEXING
            ## << TO DO >>: pre-allocate ocn_freezing_temperature
            for iCell in range(self.nCells-1):
                    kmin = self.minLevelCell[iCell]
                    kmax = self.maxLevelCell[iCell]
                    self.newPr[iTime,kmin,iCell] = self.newAtmPr[iTime,iCell] + self.newDens[iTime,iCell,kmin]*gravity*0.5*self.newLThick[iTime,iCell,kmin]

                    for k in np.arange(kmin, kmax):
                        self.newPr[iTime,k,iCell] = self.newPr[iTime,iCell,k] + 0.5*gravity*(self.newDens[iTime,iCell,k]*self.newLThick[iTime,iCell,k] \
                                + self.newDens[iTime,iCell,k-1]*self.newLThick[iTime,iCell,k-1])

                    if (self.newFloatIceMask[iTime,iCell] == 1):
                        ocn_freezing_temperature[iTime,iCell,:] = self.coeff_0_openOcean + self.coeff_S_openOcean * self.newSal[iTime,iCell,:] + self.coeff_p_openOcean * self.newPr[iTime,iCell,:] \
                                + self.coeff_pS_openOcean * self.newPr[iTime,iCell,:] * self.newSal[iTime,iCell,:] + self.coeff_mushy_openOcean * self.newSal[iTime,iCell,:] / (1.0 - self.newSal[iTime,iCell,:] / 1e3)
                    
                    elif (self.newFloatingIceMask[iTime,iCell] == 0):
                        ocn_freezing_temperature[iTime,iCell,:] = self.coeff_0_cavity + self.coeff_S_cavity * self.newSal[iTime,iCell,:] + self.coeff_p_cavity * self.newPr[iTime,iCell,:] \
                                + self.coeff_pS_cavity * self.newPr[iTime,iCell,:] * self.newSal[iTime,iCell,:] + self.coeff_mushy_cavity * self.newSal[iTime,iCell,:] / (1.0 - self.newSal[iTime,iCell,:] / 1e3)

        # Calculate thermal forcing
        self.oceanThermalForcing = self.newTemp - ocn_freezing_temperature
        
    def remap_mpas_to_mali(self)

        # make copy of mpasMeshFile to save interpolated variable to
        tmp_mpasFile = "tmp_" + self.options.mpasMeshFile
        shutil.copy(self.options.mpasMeshFile, tmp_mpasMeshFile)
        # << TO DO >>: populate tmp_mpasMeshFile with variables to be interpolated

        #create scrip files
        print("Creating Scrip Files")
        mali_scripfile = "tmp_mali_scrip.nc"
        mpas_scripfile = "tmp_mpas_scrip.nc"

        scrip_from_mpas(self.options.maliFile, mali_scripfile)
        scrip_from_mpas(tmp_mpasMeshFile, mpas_scripfile)

        #create mapping file
        print("Creating Mapping File")
        args_esmf = ['srun', 
                '-n', self.options.ntasks, 'ESMF_RegridWeightGen',
                '--source', mpas_scripfile,
                '--destination', mali_scripfile,
                '--weight', "tmp_mapping_file.nc",
                '--method', self.options.method,
                '-i', '-64bit_offset',
                "--dst_regional", "--src_regional", '--ignore_unmapped']
        check_call(args_esmf)

        #<< TO DO >>: Create mask indentifying valid overlapping ocean cells

        subprocess.run(["ncpdq", "-O", "-a", "Time,nVertLevels,nCells", tmp_mpasMeshFile, tmp_mpasMeshFile])
        subprocess.run(["ncpdq", "-O", "-a", "Time,nVertLevelsP1,nCells", tmp_mpasMeshFile, tmp_mpasMeshFile])
        subprocess.run(["ncpdq", "-O", "-a", "Time,nVertInterfaces,nCells", tmp_mpasMeshFile, tmp_mpasMeshFile])

        # remap the input data
        args_remap = ["ncremap",
                "-i", tmp_mpasMeshFile,
                "-o", self.options.outputFile,
                "-m", "tmp_mapping_file.nc"]
        check_call(args_remap)

        #remove temporary files
        os.remove("tmp_*.nc")

def main():
        run = mpasToMaliInterp()
        
        #compute yearly time average if necessary 
        run.time_average_output(self)
  
        #calculate thermal forcing
        run.calc_ocean_thermal_forcing(self)

        # << TO DO >>: copy landIceFreshwaterFlux to mpas mesh file copy for interpolation. Make name relevant to MALI

        # << TO DO >>: Vertically interpolate between MPAS-O to MALI grids

        #remap mpas to mali
        run.remap_mpas_to_mali(self)

if __name__ == "__main__":
    main()



