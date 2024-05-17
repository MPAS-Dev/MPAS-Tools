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
        print("Gathering Information ... ")
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
        IM = xr.open_dataset(self.options.maliFile, decode_times=False, decode_cf=False)
        
        # variables for time averaging
        self.temperature = DS['timeMonthly_avg_activeTracers_temperature'][:,:,:].compute()
        self.salinity = DS['timeMonthly_avg_activeTracers_temperature'][:,:,:].compute()
        self.density = DS['timeMonthly_avg_density'][:,:,:].compute()
        self.atmPressure = DS['timeMonthly_avg_atmosphericPressure'][:,:].compute()
        self.daysSinceStart = DS['timeMonthly_avg_daysSinceStartOfSim'][:].compute()
        self.stTime = DS.attrs['config_start_time']
        self.landIceFloatingMask = OM['landIceFloatingMask'][0,:].data #not letting floating ice mask evolve for now because it's only in the mesh file
        self.layerThickness = DS['timeMonthly_avg_layerThickness'][:,:].data
        xt = DS['xtime_startMonthly']
        xtime = np.array([xt],dtype=('S',64)) 
        # << NOTE >>: may need to eventually use: "xtime.data.tobytes().decode()" but not sure yet
        try:
            self.landIceFreshwaterFlux = DS['timeMonthly_avg_landIceFreshwaterFlux'][:,:].data
        except KeyError:
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
        az1_liq = DS.attrs['config_open_ocean_freezing_temperature_coeff_mushy_az1_liq']
        self.coeff_mushy_openOcean = 1/az1_liq 
        self.coeff_0_cavity = DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_0']
        self.coeff_S_cavity = DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_S']
        self.coeff_p_cavity = DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_p']
        self.coeff_pS_cavity = DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_pS']
        self.coeff_mushy_cavity = 0

        #variables for interpolation << NOTE >>: a couple of these fields are redundant but we want them in xarray form for compute_zmid
        self.bottomDepth_dataArray = OM['bottomDepth'][:]
        self.maxLevelCell_dataArray = OM['maxLevelCell'][:]
        self.layerThickness_dataArray = DS['timeMonthly_avg_layerThickness']
        
    def time_average_output(self):
        print("Time Averaging ...")
        if (self.options.yearlyAvg == 'true'):

            print("daysSinceStart: {}".format(self.daysSinceStart))
            yearsSinceStart = self.daysSinceStart / 365.0
            print("yearsSinceStart: {}".format(yearsSinceStart))
            finalYear = np.floor(np.max(yearsSinceStart))
            startYear = np.floor(np.min(yearsSinceStart))
            print("final year: {}".format(finalYear))
            print("start year: {}".format(startYear))
            timeStride = 1 # 1 year average
            
            if (startYear != finalYear):
                years = np.arange(startYear, finalYear, timeStride, dtype=int)
                nt = len(years)
            else :
                years = startYear
                nt = 1

            print("years {}: ".format(years))

            #pre-allocate
            _,nc,nz = self.temperature.shape
            self.newTemp = np.zeros((nt,nc,nz)) 
            self.newSal = np.zeros((nt,nc,nz))
            self.newDens = np.zeros((nt,nc,nz))
            self.newLThick= np.zeros((nt,nc,nz))
            self.newAtmPr = np.zeros((nt,nc))
            if (self.landIceFreshwaterFlux != 0):
                self.newLandIceFWFlux = np.zeros((nt,nc))

            print("starting loop ...")
            if (years.ndim == 0): 
                print("Years: {}".format(years))
                log = np.logical_and(yearsSinceStart >= years, yearsSinceStart < years + timeStride)
                self.newTemp[0,:,:] = np.mean(self.temperature[log,:,:], axis=0)
                self.newSal[0,:,:] = np.mean(self.salinity[log,:,:], axis=0)
                self.newDens[0,:,:] = np.mean(self.density[log,:,:], axis=0)
                self.newLThick[0,:,:] = np.mean(self.layerThickness[log,:,:], axis=0)
                self.newAtmPr[0,:] = np.mean(self.atmPressure[log,:], axis=0)
                if (self.landIceFreshwaterFlux != 0):
                    self.newLandIceFWFlux[0,:] = np.mean(self.landIceFreshwaterFlux[log,:], axis=0)
                else :
                    self.newLandIceFWFlux = 0
            else :
                ct = 0
                for i in years:
                    print("Year: {}".format(years[i]))
                    log = np.logical_and(yearsSinceStart >= years[i], yearsSinceStart < years[i] + timeStride)
                    self.newTemp[ct,:,:] = np.mean(self.temperature[log,:,:], axis=0)
                    self.newSal[ct,:,:] = np.mean(self.salinity[log,:,:], axis=0)
                    self.newDens[ct,:,:] = np.mean(self.density[log,:,:], axis=0)
                    self.newLThick[ct,:,:] = np.mean(self.layerThickness[log,:,:], axis=0)
                    self.newAtmPr[ct,:] = np.mean(self.atmPressure[log,:], axis=0)
                    if (self.landIceFreshwaterFlux != 0):
                        self.newLandIceFWFlux[ct,:] = np.mean(self.landIceFreshwaterFlux[log,:], axis=0)
                    else :
                        self.newLandIceFWFlux = 0
                    ct = ct + 1

            # Build new xtime string. Defined on first day of each year
            stTime = np.datetime64(self.stTime[0:4])
            stYear = np.datetime_as_string(stTime,unit='Y').astype('float64')
            yearVec = yearsSinceStart + stYear
            self.newXtime = np.array([f"{int(year)}-01-01_00:00:00" for year in yearVec],dtype=np.str_)
        else : #do nothing if not time averaging
            self.newTemp = self.temperature
            self.newSal = self.salinity
            self.newDens= self.density
            self.newLThick = self.layerThickness
            self.newAtmPr = self.atmPressure
            self.newLandIceFWFlux = self.landIceFreshwaterFlux
            self.newXtime = xtime
    
    def calc_ocean_thermal_forcing(self):
        print("Calculating thermal forcing ... ")
        gravity = 9.81

        #pre-allocate
        nt,nc,nz = self.newTemp.shape
        ocn_freezing_temperature = np.zeros((nt,nc,nz))
        self.newPr = np.zeros((nt,nc,nz))

        print("floating ice mask: {}".format(self.landIceFloatingMask.shape))
        for iTime in range(nt):

            #calculate pressure: 
            for iCell in range(self.nCells):
                    kmin = self.minLevelCell[iCell] - 1
                    kmax = self.maxLevelCell[iCell] - 1
                    self.newPr[iTime,kmin,iCell] = self.newAtmPr[iTime,iCell] + self.newDens[iTime,iCell,kmin]*gravity*0.5*self.newLThick[iTime,iCell,kmin]

                    for k in np.arange(kmin + 1, kmax):
                        self.newPr[iTime,k,iCell] = self.newPr[iTime,iCell,k-1] + 0.5*gravity*(self.newDens[iTime,iCell,k-1]*self.newLThick[iTime,iCell,k-1] \
                                + self.newDens[iTime,iCell,k]*self.newLThick[iTime,iCell,k])

                    print("floating ice mask iCell: {}".format(self.landIceFloatingMask[iCell]))
                    if (self.landIceFloatingMask[iCell] == 1):
                        ocn_freezing_temperature[iTime,iCell,:] = self.coeff_0_openOcean + self.coeff_S_openOcean * self.newSal[iTime,iCell,:] + self.coeff_p_openOcean * self.newPr[iTime,iCell,:] \
                                + self.coeff_pS_openOcean * self.newPr[iTime,iCell,:] * self.newSal[iTime,iCell,:] + self.coeff_mushy_openOcean * self.newSal[iTime,iCell,:] / (1.0 - self.newSal[iTime,iCell,:] / 1e3)
                    
                    elif (self.landIceFloatingMask[iCell] == 0):
                        ocn_freezing_temperature[iTime,iCell,:] = self.coeff_0_cavity + self.coeff_S_cavity * self.newSal[iTime,iCell,:] + self.coeff_p_cavity * self.newPr[iTime,iCell,:] \
                                + self.coeff_pS_cavity * self.newPr[iTime,iCell,:] * self.newSal[iTime,iCell,:] + self.coeff_mushy_cavity * self.newSal[iTime,iCell,:] / (1.0 - self.newSal[iTime,iCell,:] / 1e3)

        # Calculate thermal forcing
        self.oceanThermalForcing = self.newTemp - ocn_freezing_temperature
        
    def remap_mpas_to_mali(self):
        print("Start remapping ... ")

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

        print("ncremap ...")
        # remap the input data
        args_remap = ["ncremap",
                "-i", tmp_mpasMeshFile,
                "-o", "tmp_" + self.options.outputFile,
                "-m", "tmp_mapping_file.nc"]
        check_call(args_remap)

        print("Finishing processing and saving ...")
        # Create mask indentifying valid overlapping ocean cells. Combine all MALI terms in one input file
        interpDS = xr.open_dataset("tmp_" + self.options.outputFile, decode_times=False, decode_cf=False)
        interpTF = interpDS['ismip6shelfMelt_3dThermalForcing'][:,:,:]
        if (self.newLandIceFWFlux != 0):
            interpFWF = interpDS['floatingBasalMassBal'][:,:]

        validOceanMask = np.zeros(self.interpTF.data.shape, dtype='int')
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
        run.time_average_output()
  
        #calculate thermal forcing
        run.calc_ocean_thermal_forcing()

        #remap mpas to mali
        run.remap_mpas_to_mali()

if __name__ == "__main__":
    main()



