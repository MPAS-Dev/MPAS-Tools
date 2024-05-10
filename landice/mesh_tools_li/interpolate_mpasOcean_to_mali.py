#!/usr/bin/env python

import os
import shutil
import xarray as xr
import subprocess
import numpy as np
from scipy import stats as st
from mpas_tools.logging import check_call
from mpas_tools.scrip.from_mpas import scrip_from_mpas
from mpas_tools.ocean.depth import compute_zmid
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
        options, args = parser.parse_args()
        self.options = options

        # open and concatenate diagnostic dataset
        DS = xr.open_mfdataset(self.options.mpasDiagsDir + '/' + '*.nc', combine='nested', concat_dim='Time', decode_timedelta=False)        
        #open mpas ocean mesh
        OM = xr.open_dataset(self.options.mpasMeshFile, decode_times=False, decode_cf=False)
        #open MALI mesh
        IM = xr.open_dataset(self.options.maliMeshFile, decode_times=False, decode_cf=False)
        
        # variables for time averaging
        self.temperature = DS['timeMonthly_avg_activeTracers_temperature'][:,:,:].data
        self.salinity = DS['timeMonthly_avg_activeTracers_temperature'][:,:,:].data
        self.density = DS['timeMonthly_avg_density'][:,:,:].data
        self.daysSinceStart = DS['timeMonthly_avg_daysSinceStartOfSim'][:].data
        self.stTime = DS.attrs['config_start_time']
        self.landIceFloatingMask = OM['landIceFloatingMask'][:,:,:].data
        self.layerThickness = DS['timeMonthly_avg_layerThickness'][:,:].data
        xt = DS['xtime'][:,:].data
        xtime = np.array([xt],dtype=('S',64)) 
        # << NOTE >>: may need to eventually use: "xtime.data.tobytes().decode()" but not sure yet
        try:
            self.landIceFreshwaterFlux = DS['timeMonthly_avg_landIceFreshwaterFlux'][:,:].data
        finally:
            print("No landIceFreshwaterFlux variable")
            self.landIceFreshwaterFlux = 0

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

        #variables for interpolation << NOTE >>: a couple of these fields are redundant but we want them in xarray form for compute_zmid
        self.bottomDepth_dataArray = OM['bottomDepth'][:]
        self.maxLevelCell_dataArray = OM['maxLevelCell'][:]
        self.layerThickness_dataArray = DS['timeMonthly_avg_layerThickness']
        
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
                if (self.landIceFreshwaterFlux != 0):
                    self.newLandIceFWFlux[i,:] = np.mean(self.landIceFreshwaterFlux[ind,:], axis=0)
                elif
                    self.newLandIceFWFlux = 0

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
            self.newLandIceFWFlux = self.landIceFreshwaterFlux
            self.newXtime = xtime
    
    def calc_ocean_thermal_forcing(self)

        gravity = 9.81
        nt,nc,nz = self.newTemp.shape

        ocn_freezing_temperature = np.zeros((nt,nc,nz))
        for iTime in range(nt):

            #calculate pressure: 
            for iCell in range(self.nCells):
                    kmin = self.minLevelCell[iCell] - 1
                    kmax = self.maxLevelCell[iCell] - 1
                    self.newPr[iTime,kmin,iCell] = self.newAtmPr[iTime,iCell] + self.newDens[iTime,iCell,kmin]*gravity*0.5*self.newLThick[iTime,iCell,kmin]

                    for k in np.arange(kmin + 1, kmax):
                        self.newPr[iTime,k,iCell] = self.newPr[iTime,iCell,k-1] + 0.5*gravity*(self.newDens[iTime,iCell,k-1]*self.newLThick[iTime,iCell,k-1] \
                                + self.newDens[iTime,iCell,k]*self.newLThick[iTime,iCell,k])

                    if (self.newFloatingIceMask[iTime,iCell] == 1):
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
        
        # populate tmp_mpasMeshFile with variables to be interpolated
        f = xr.open_dataset(tmp_mpasMeshFile, decode_times=False, decode_cf=False)
        tf = xr.DataArray(self.oceanThermalForcing.astype('float64'),dims=('Time','nCells','nVertLevels'))
        f['ismip6shelfMelt_3dThermalForcing'] = tf

        if (self.newLandIceFWFlux != 0):
            m = xr.DataArray(self.newLandIceFWFlux('float64'),dims=('Time','nCells'))
            f['floatingBasalMassBal'] = m

        # save new mesh file
        f.to_netcdf(tmp_mpasMeshFile)
        f.close()

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

        subprocess.run(["ncpdq", "-O", "-a", "Time,nVertLevels,nCells", tmp_mpasMeshFile, tmp_mpasMeshFile])
        subprocess.run(["ncpdq", "-O", "-a", "Time,nVertLevelsP1,nCells", tmp_mpasMeshFile, tmp_mpasMeshFile])
        subprocess.run(["ncpdq", "-O", "-a", "Time,nVertInterfaces,nCells", tmp_mpasMeshFile, tmp_mpasMeshFile])

        # remap the input data
        args_remap = ["ncremap",
                "-i", tmp_mpasMeshFile,
                "-o", "tmp_" + self.options.outputFile,
                "-m", "tmp_mapping_file.nc"]
        check_call(args_remap)

        # Create mask indentifying valid overlapping ocean cells. Combine all MALI terms in one input file
        interpDS = xr.open_dataset("tmp_" + self.options.outputFile, decode_times=False, decode_cf=False)
        interpTF = interpDS['ismip6shelfMelt_3dThermalForcing'][:,:,:]
        if (self.newLandIceFWFlux != 0):
            interpFWF = interpDS['floatingBasalMassBal'][:,:]

        validOceanMask = np.zeros(self.interpTF.data.shape, np.dtype = 'int')
        ind = np.where(interpTF != 0)
        validOceanMask[ind] = 1

        mask = xr.dataArray(validOceanMask,dims=("Time","nCells"),long_name="Mask of MALI cells overlapping interpolated MPAS-Oceans cells")
        IM = open_dataset(self.options.maliFile,decode_times=False,decode_cf=False)
        IM['validOceanMask'] = mask
        if (self.newLandIceFWFlux != 0):
            IM['floatingBasalMassBal'] = interpFWF
        IM['ismip6shelfMelt_3dThermalForcing'] = interpTF

        # define zmid coordinate in MPAS-O and translate to ISMIP6 vertical coordinate
        zmid = compute_zmid(self.bottomDepth_dataArray, self.maxLevelCell_dataArray, self.layerThickness_dataArray)
        IM['ismip6shelfMelt_zOcean'] = zmid

        IM.to_netcdf(self.options.outputFile)
        IM.close()

        #remove temporary files
        os.remove("tmp_*.nc")

def main():
        run = mpasToMaliInterp()
        
        #compute yearly time average if necessary 
        run.time_average_output(self)
  
        #calculate thermal forcing
        run.calc_ocean_thermal_forcing(self)

        #remap mpas to mali
        run.remap_mpas_to_mali(self)

if __name__ == "__main__":
    main()



