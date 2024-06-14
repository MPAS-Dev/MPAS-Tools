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
import time

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
        avg_layerThickness = DS['timeMonthly_avg_layerThickness']
        self.layerThickness = avg_layerThickness.data
        alT = avg_layerThickness.compute().values
        print("avg_layerThickness: {}".format(alT.shape))
        print("avg_layerTHickness NaNs: {}".format(np.sum(np.isnan(alT))))

        xt = DS['xtime_startMonthly']
        xtime = np.array([xt],dtype=('S',64)) 
        # << NOTE >>: may need to eventually use: "xtime.data.tobytes().decode()" but not sure yet
        try:
            self.landIceFreshwaterFlux = DS['timeMonthly_avg_landIceFreshwaterFlux'][:,:].data
        except KeyError:
            print("No landIceFreshwaterFlux variable")
            self.landIceFreshwaterFlux = 0

        #check for nans in original data
        if np.isnan(self.temperature).any():
            print("Original temperature contains NaNs")
        else:
            print("Original temperature does not contain NaNs")

        # additional variables for computing thermal forcing 
        bottomDepth = OM['bottomDepth']
        bD = bottomDepth.data
        print("bottomDepth: {}".format(bD.shape))
        print("bottomDepth NaNs: {}".format(np.sum(np.isnan(bD))))

        self.minLevelCell = OM['minLevelCell'][:].data
        maxLevelCell = OM['maxLevelCell']
        self.maxLevelCell = maxLevelCell.data
        print("maxLevelCell: {}".format(self.maxLevelCell.shape))
        print("maxLevelCell NaNs: {}".format(np.sum(np.isnan(self.maxLevelCell))))

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

        # Define vertical coordinates of mpas output
        mpas_cellCenterElev = compute_zmid(bottomDepth, maxLevelCell, avg_layerThickness)
        self.mpas_cellCenterElev = mpas_cellCenterElev.data
        print("mpas_cellCenterElev {}".format(self.mpas_cellCenterElev.shape))
        mpas_cellCenterElev = mpas_cellCenterElev.compute().values
        print("mpas_cellCenterElev number of NaNs: {}".format(np.sum(np.isnan(mpas_cellCenterElev))))

    def time_average_output(self):
        print("Time Averaging ...")
        if (self.options.yearlyAvg == 'true'):

            yearsSinceStart = self.daysSinceStart / 365.0
            finalYear = np.floor(np.max(yearsSinceStart))
            startYear = np.floor(np.min(yearsSinceStart))
            timeStride = 1 # 1 year average
            
            if (startYear != finalYear):
                years = np.arange(startYear, finalYear, timeStride, dtype=int)
                nt = len(years)
            else :
                years = startYear
                nt = 1

            #pre-allocate
            _,nc,nz = self.temperature.shape
            self.newTemp = np.zeros((nt,nc,nz)) 
            self.newSal = np.zeros((nt,nc,nz))
            self.newDens = np.zeros((nt,nc,nz))
            self.newLThick= np.zeros((nt,nc,nz))
            self.newAtmPr = np.zeros((nt,nc))
            self.newMpasCCE = np.zeros((nt,nc,nz))
            if (self.landIceFreshwaterFlux != 0):
                self.newLandIceFWFlux = np.zeros((nt,nc))

            print("starting loop ...")
            st = time.time()
            if (years.ndim == 0): 
                log = np.logical_and(yearsSinceStart >= years, yearsSinceStart < years + timeStride)
                self.newTemp[0,:,:] = np.mean(self.temperature[log,:,:], axis=0)
                self.newSal[0,:,:] = np.mean(self.salinity[log,:,:], axis=0)
                self.newDens[0,:,:] = np.mean(self.density[log,:,:], axis=0)
                self.newLThick[0,:,:] = np.mean(self.layerThickness[log,:,:], axis=0)
                self.newMpasCCE[0,:,:] = np.mean(self.mpas_cellCenterElev[log,:,:], axis=0)
                self.newAtmPr[0,:] = np.mean(self.atmPressure[log,:], axis=0)
                if (self.landIceFreshwaterFlux != 0):
                    self.newLandIceFWFlux[0,:] = np.mean(self.landIceFreshwaterFlux[log,:], axis=0)
                else :
                    self.newLandIceFWFlux = 0
            else :
                ct = 0
                for i in years:
                    log = np.logical_and(yearsSinceStart >= years[i], yearsSinceStart < years[i] + timeStride)
                    self.newTemp[ct,:,:] = np.mean(self.temperature[log,:,:], axis=0)
                    self.newSal[ct,:,:] = np.mean(self.salinity[log,:,:], axis=0)
                    self.newDens[ct,:,:] = np.mean(self.density[log,:,:], axis=0)
                    self.newLThick[ct,:,:] = np.mean(self.layerThickness[log,:,:], axis=0)
                    self.newMpasCCE[ct,:,:] = np.mean(self.mpas_cellCenterElev[log,:,:], axis=0)
                    self.newAtmPr[ct,:] = np.mean(self.atmPressure[log,:], axis=0)
                    if (self.landIceFreshwaterFlux != 0):
                        self.newLandIceFWFlux[ct,:] = np.mean(self.landIceFreshwaterFlux[log,:], axis=0)
                    else :
                        self.newLandIceFWFlux = 0
                    ct = ct + 1
            nd = time.time()
            tm = nd - st
            print("Time averaging loop", tm, "seconds")

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
            self.newMpasCCE = self.mpas_cellCenterElev
            self.newAtmPr = self.atmPressure
            self.newLandIceFWFlux = self.landIceFreshwaterFlux
            self.newXtime = xtime

        print("neMPASCCE: {}".format(self.newMpasCCE.shape))
        print("NaNs in newMpasCCE: {}".format(np.sum(np.isnan(self.newMpasCCE))))

        #check for nans
        if np.isnan(self.mpas_cellCenterElev).any():
            print("mpas_cellCenterElev contains NaNs")
        else:
            print("mpas_cellCenterElev does not contain NaNs")

        if np.isnan(self.newMpasCCE).any():
            print("newMpasCCE contains NaNs")
        else:
            print("newMpasCCE does not contain NaNs")
        
        if np.isnan(self.temperature).any():
            print("temperature contains NaNs")
        else:
            print("temperature does not contain NaNs")
    
        if np.isnan(self.newTemp).any():
            print("newTemp contains NaNs")
        else:
            print("newTemp does not contain NaNs")
    
    def calc_ocean_thermal_forcing(self):
        print("Calculating thermal forcing ... ")
        gravity = 9.81

        #pre-allocate
        nt,nc,nz = self.newTemp.shape
        ocn_freezing_temperature = np.zeros((nt,nc,nz))
        self.newPr = np.zeros((nt,nc,nz))
        st = time.time()
        for iTime in range(nt):

            #calculate pressure: 
            for iCell in range(self.nCells):
                    kmin = self.minLevelCell[iCell] - 1
                    kmax = self.maxLevelCell[iCell] - 1

                    self.newPr[iTime,iCell,kmin] = self.newAtmPr[iTime,iCell] + \
                            self.newDens[iTime,iCell,kmin]*gravity*0.5*self.newLThick[iTime,iCell,kmin]

                    for k in np.arange(kmin + 1, kmax):
                        self.newPr[iTime,iCell,k] = self.newPr[iTime,iCell,k-1] + \
                                0.5*gravity*(self.newDens[iTime,iCell,k-1]*self.newLThick[iTime,iCell,k-1] \
                                + self.newDens[iTime,iCell,k]*self.newLThick[iTime,iCell,k])

                    if (self.landIceFloatingMask[iCell] == 1):
                        ocn_freezing_temperature[iTime,iCell,:] = self.coeff_0_openOcean + \
                                self.coeff_S_openOcean * self.newSal[iTime,iCell,:] + self.coeff_p_openOcean * self.newPr[iTime,iCell,:] \
                                + self.coeff_pS_openOcean * self.newPr[iTime,iCell,:] * self.newSal[iTime,iCell,:] \
                                + self.coeff_mushy_openOcean * self.newSal[iTime,iCell,:] / (1.0 - self.newSal[iTime,iCell,:] / 1e3)
                    
                    elif (self.landIceFloatingMask[iCell] == 0):
                        ocn_freezing_temperature[iTime,iCell,:] = self.coeff_0_cavity + \
                                self.coeff_S_cavity * self.newSal[iTime,iCell,:] + self.coeff_p_cavity * self.newPr[iTime,iCell,:] \
                                + self.coeff_pS_cavity * self.newPr[iTime,iCell,:] * self.newSal[iTime,iCell,:] \
                                + self.coeff_mushy_cavity * self.newSal[iTime,iCell,:] / (1.0 - self.newSal[iTime,iCell,:] / 1e3)
        nd = time.time()
        tm = nd - st
        print("Ocean thermal forcing loop", tm, "seconds")
        # Calculate thermal forcing
        self.oceanThermalForcing = self.newTemp - ocn_freezing_temperature
        
        if np.isnan(self.oceanThermalForcing).any():
            print("oceanThermalForcing contains NaNs")
        else:
            print("oceanThermalForcing does not contain NaNs")
    
    def remap_mpas_to_mali(self):
        print("Start remapping ... ")

        # populate tmp_mpasMeshFile with variables to be interpolated
        f = xr.open_dataset(self.options.mpasMeshFile, decode_times=False, decode_cf=False)
        tf = xr.DataArray(self.oceanThermalForcing.astype('float64'),dims=('Time','nCells','nVertLevels'))
        
        f['ismip6shelfMelt_3dThermalForcing'] = tf

        mcce = xr.DataArray(self.newMpasCCE.astype('float64'),dims=('Time','nCells','nVertLevels'))
        f['mpas_cellCenterElev'] = mcce

        if (self.newLandIceFWFlux != 0):
            m = xr.DataArray(self.newLandIceFWFlux.astype('float64'),dims=('Time','nCells'))
            f['floatingBasalMassBal'] = m
        
        #delete unncessary variables with 'string1' dimension or will cause errors.
        try:
            f = f.drop_vars('simulationStartTime')
        except ValueError:
            f = f 

        try:
            f = f.drop_vars('forcingGroupNames')
        except ValueError:
            f = f 
        
        try:
            f = f.drop_vars('forcingGroupRestartTimes')
        except ValueError:
            f = f 
        
        st = time.time()
        # save new mesh file
        tmp_mpasMeshFile = "tmp_mpasMeshFile.nc"
        f.to_netcdf(tmp_mpasMeshFile)
        f.close()
        nd = time.time()
        tm = nd - st
        print("saving updates mpas mesh file", tm, "seconds")

        st = time.time()
        subprocess.run(["ncatted", "-a", "_FillValue,,d,,", tmp_mpasMeshFile])
        nd = time.time()
        tm = nd - st
        print("Removing fill value:", tm, "seconds")

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
                "-o", "tmp_outputFile.nc",
                "-m", "tmp_mapping_file.nc"]
        check_call(args_remap)

        #Vertical interpolation
        print("Vertically interpolating onto mali grid")
        interpDS = xr.open_dataset("tmp_outputFile.nc", decode_times=False, decode_cf=False)
        TF = interpDS['ismip6shelfMelt_3dThermalForcing'][:,:,:].data
        mpas_cellCenterElev = interpDS['mpas_cellCenterElev'][:,:,:].data

        #Define vertical coordinates for mali grid
        IM = xr.open_dataset(self.options.maliFile,decode_times=False,decode_cf=False)
        layerThicknessFractions = IM['layerThicknessFractions'][:].data
        bedTopography = IM['bedTopography'][:,:].data
        thickness = IM['thickness'][:,:].data
        
        nz, = layerThicknessFractions.shape
        nt,nc = thickness.shape

        layerThicknessFractions = layerThicknessFractions.reshape((1,1,nz))
        thickness = thickness.reshape((nt,nc,1))
        bedTopography = bedTopography.reshape((nt,nc,1))

        layerThickness = thickness * layerThicknessFractions
        mali_cellCenterElev = np.cumsum(layerThickness,axis=2) - layerThickness/2 + bedTopography

        # Reshape to (nt*nc, nz) before interpolation to avoid looping
        mali_cellCenterElev = mali_cellCenterElev.reshape(nt*nc, -1)
        TF = np.transpose(TF,(0,2,1))
        TF = TF.reshape(nt*nc, -1)
        mpas_cellCenterElev = np.transpose(mpas_cellCenterElev,(0,2,1))
        mpas_cellCenterElev = mpas_cellCenterElev.reshape(nt*nc, -1)

        print("Before Interpolation")
        print("TF zeroth dim: {}".format(TF[0:5,0]))
        print("TF first dim: {}".format(TF[0,0:5]))

        print("mpas_cellCenterElev zeroth dim: {}".format(mpas_cellCenterElev[0:5,0]))
        print("mpas_cellCenterElev first dim: {}".format(mpas_cellCenterElev[0,0:5]))

        st = time.time()
        # Linear interpolation 
        
        if np.isnan(mali_cellCenterElev).any():
            print("mali_cellCenterElev contains NaNs")
        else:
            print("mali_cellCenterElev does not contain NaNs")

        if np.isnan(mpas_cellCenterElev).any():
            print("mpas_cellCenterElev contains NaNs")
        else:
            print("mpas_cellCenterElev does not contain NaNs")

        if np.isnan(TF).any():
            print("TF contains NaNs")
        else:
            print("TF does not contain NaNs")

        nan_count = 0
        for i in range(nt*nc):
            ind = np.where(~np.isnan(TF[i,:].flatten() * mpas_cellCenterElev[i,:].flatten()))
            
            #if (ind.size != 0): #only interpolation where valid data
            ct = 0
            if len(ind[0]) != 0:
                if (ct == 0):
                    print("mali_cellCenterElev[i,:]".format(mali_cellCenterElev[i,:].flatten().shape))
                    print("mpas_cellCenterElev[i,ind]".format(mpas_cellCenterElev[i,ind].flatten().shape))
                    print("TF[i,ind]".format(mali_cellCenterElev[i,ind].flatten().shape))
                    ct = ct + 1
                interpTF = np.array(np.interp(mali_cellCenterElev[i,:].flatten(), mpas_cellCenterElev[i,ind].flatten(), TF[i,ind].flatten()))
            else:
                nan_count = nan_count + 1
        print("NaN Count: {}; nt*nc = {}".format(nan_count, nt*nc))
        nd = time.time()
        tm = nd-st
        print("vertical interpolation time: {}",tm, "seconds")

        print("After interpolation:")
        print("interpTF: {}".format(interpTF[0:5,0]))
        
        # Reshape interpTF back to (nt, nc, nz)
        interpTF = interpTF.reshape(nt, nc, -1)

        print("After Reshape:")
        print("interpTF: {}".format(interpTF[0,0:5,0]))

        # Create mask indentifying valid overlapping ocean cells. Combine all MALI terms in one input file
        #Making this a 3-D variable for now, but may want to make 2-D eventually
        validOpenOceanMask = np.zeros((nt,nc), dtype='int')
        ind = np.where(interpTF[:,:,0].data != 0)
        validOpenOceanMask[ind] = 1

        if (self.newLandIceFWFlux != 0):
            interpFWF = interpDS['floatingBasalMassBal'][:,:]
            print("interpFWF: {}".format(interpFWF.data.shape))

        print("Saving")
        
        print("interpTF: {}".format(interpTF.shape))
        print("validOpenOceanMask: {}".format(validOpenOceanMask.shape))

        mask = xr.DataArray(validOpenOceanMask,dims=("Time","nCells"),attrs={'long_name':'Mask of MALI cells overlapping interpolated MPAS-Oceans cells'})
        interptf = xr.DataArray(interpTF, dims=("Time","nCells","nVertLevels"))

        IM['validOceanMask'] = mask
        if (self.newLandIceFWFlux != 0):
            IM['floatingBasalMassBal'] = interpFWF
        IM['ismip6shelfMelt_3dThermalForcing'] = interptf

        IM.to_netcdf(self.options.outputFile)
        IM.close()

#        subprocess.run(["ncatted", "-a", "_FillValue,,d,,", self.options.outputFile])

        #remove temporary files
        files = os.listdir()
        for file in files:
            if file.startswith("tmp"):
                os.remove(file)

def main():
        run = mpasToMaliInterp()
        
        #compute yearly time average if necessary 
        run.time_average_output()
  
        #calculate thermal forcing
        run.calc_ocean_thermal_forcing()

        #remap mpas to mali
        run.remap_mpas_to_mali()
        
        print("Finished.")
if __name__ == "__main__":
    main()



