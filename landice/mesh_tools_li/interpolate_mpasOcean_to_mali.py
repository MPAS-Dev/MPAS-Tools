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
import cftime

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

        #xt = DS['xtime_startMonthly']
        
        #xtime = np.array([xt],dtype=('S',64)) 
        # << NOTE >>: may need to eventually use: "xtime.data.tobytes().decode()" but not sure yet
        
        try:
            self.landIceFreshwaterFlux = DS['timeMonthly_avg_landIceFreshwaterFlux'][:,:].data
        except KeyError:
            print("No landIceFreshwaterFlux variable")
            self.landIceFreshwaterFlux = 0

        # additional variables for computing thermal forcing 
        bottomDepth = OM['bottomDepth']
        bD = bottomDepth.data

        self.minLevelCell = OM['minLevelCell'][:].data
        maxLevelCell = OM['maxLevelCell']
        self.maxLevelCell = maxLevelCell.data

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

    def time_average_output(self):
        print("Time Averaging ...")
        if (self.options.yearlyAvg == 'true'):

            yearsSinceStart = self.daysSinceStart.data / 365.0
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
            yearVec = np.zeros((nt,))    
            #prepare time vector
            stTime = np.datetime64(self.stTime[0:4])
            print("stTime = {}".format(stTime))
            stYear = np.datetime_as_string(stTime,unit='Y').astype('float64')
            print("stYear = {}".format(stYear))
            
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
            
                #Define time at the first of each year

                print("Test A: {}".format(np.floor(np.min(yearsSinceStart[log])) + stYear))
                yearVec[0] = np.floor(np.min(yearsSinceStart[log])) + stYear
                print("yearVec = {}".format(yearVec))
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
       
                    yearVec[ct] = np.floor(np.min(yearsSinceStart[log])) + stYear
                    print("yearVec = {}".format(yearVec))
                    ct = ct + 1
            nd = time.time()
            tm = nd - st
            print("Time averaging loop", tm, "seconds")
            
            #establish xtime
            dates = []
            xt = cftime.num2date(yearVec*365, units="days since 0000-01-01_00:00:00",calendar='noleap')
            xt_reformat = [dt.strftime("%Y-%m-%d_%H:%M:%S") for dt in xt]
            dates.append(xt_reformat)
            #self.newXtime = dates        
            self.newXtime = xt_reformat
        else : #do nothing if not time averaging
            self.newTemp = self.temperature
            self.newSal = self.salinity
            self.newDens= self.density
            self.newLThick = self.layerThickness
            self.newMpasCCE = self.mpas_cellCenterElev
            self.newAtmPr = self.atmPressure
            self.newLandIceFWFlux = self.landIceFreshwaterFlux
            self.newXtime = np.nan #use as placeholder for now

        print("New xtime: {}".format(self.newXtime))
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
        #For now, hardcode in depths from ISMIP6 AIS. Should work fine for our purposes
        ismip6shelfMelt_zOcean = np.array((-30, -90, -150, -210, -270, -330, -390, -450, -510,\
        -570, -630, -690, -750, -810, -870, -930, -990, -1050, -1110, -1170,\
        -1230, -1290, -1350, -1410, -1470, -1530, -1590, -1650, -1710, -1770),dtype='float64')
        
        nz, = ismip6shelfMelt_zOcean.shape

        # Reshape to (nt*nc, nz) before interpolation to avoid looping
        TF = np.transpose(TF,(0,2,1))
        mpas_cellCenterElev = np.transpose(mpas_cellCenterElev,(0,2,1))
        nt,nc,_ = mpas_cellCenterElev.shape
        print("mpas_cellCenterElev: {}".format(mpas_cellCenterElev.shape))
        print("TF: {}".format(TF.shape))
        print("nz: {}:".format(nz))
        TF = TF.reshape(nt*nc, -1)
        mpas_cellCenterElev = mpas_cellCenterElev.reshape(nt*nc, -1)

        # Linear interpolation 
        st = time.time()
        nan_count = 0
        interpTF = np.zeros((nt*nc,nz))
        for i in range(nt*nc):
            ind = np.where(~np.isnan(TF[i,:].flatten() * mpas_cellCenterElev[i,:].flatten()))
            
            ct = 0
            if len(ind[0]) != 0:
                interpTF[i,:] = np.array(np.interp(ismip6shelfMelt_zOcean, mpas_cellCenterElev[i,ind].flatten(), TF[i,ind].flatten()))
            else:
                nan_count = nan_count + 1
        nd = time.time()
        tm = nd-st
        print("vertical interpolation time: {}",tm, "seconds")

        # Reshape interpTF back to (nt, nc, nz)
        interpTF = interpTF.reshape(nt, nc, -1)

        # Create mask indentifying valid overlapping ocean cells. Combine all MALI terms in one input file
        #Making this a 3-D variable for now, but may want to make 2-D eventually
        validOpenOceanMask = np.zeros((nt,nc), dtype='int')
        print("interpTF: {}".format(interpTF.data.shape))
        surfaceTF = interpTF[:,:,0]
        ind = np.where(surfaceTF != 0)
        validOpenOceanMask[ind] = 1


        print("Saving")
        #open original mali file
        IM = xr.open_dataset(self.options.maliFile,decode_times=False,decode_cf=False)
        
        #create new mali forcing file that only desired forcing variables
        ds_out = IM.copy(deep=False)
        ds_out = ds_out.drop_vars(ds_out.data_vars)
        
        # introduce new ismip6 depth dimension
        ds_out = ds_out.expand_dims(nISMIP6OceanLayers=ismip6shelfMelt_zOcean) # introduce new ismip6 depth dimension

        # Save variables
        mask = xr.DataArray(validOpenOceanMask,dims=("Time","nCells"),attrs={'long_name':'Mask of MALI cells overlapping interpolated MPAS-Oceans cells'})       
        interptf = xr.DataArray(interpTF, dims=("Time","nCells","nISMIP6OceanLayers"))
        zlayers = xr.DataArray(ismip6shelfMelt_zOcean, dims=("nISMIP6OceanLayers"))

        ds_out['angleEdge'] = IM['angleEdge']
        ds_out['areaCell'] = IM['areaCell']
        ds_out['areaTriangle'] = IM['areaTriangle']
        ds_out['cellsOnCell'] = IM['cellsOnCell']
        ds_out['cellsOnEdge'] = IM['cellsOnEdge']
        ds_out['cellsOnVertex'] = IM['cellsOnVertex']
        ds_out['dcEdge'] = IM['dcEdge']
        ds_out['dvEdge'] = IM['dvEdge']
        ds_out['edgesOnCell'] = IM['edgesOnCell']
        ds_out['edgesOnEdge'] = IM['edgesOnEdge']
        ds_out['edgesOnVertex'] = IM['edgesOnVertex']
        ds_out['indexToCellID'] = IM['indexToCellID']
        ds_out['indexToEdgeID'] = IM['indexToEdgeID']
        ds_out['indexToVertexID'] = IM['indexToVertexID']
        ds_out['kiteAreasOnVertex'] = IM['kiteAreasOnVertex']
        ds_out['latCell'] = IM['latCell']
        ds_out['latEdge'] = IM['latEdge']
        ds_out['latVertex'] = IM['latVertex']
        ds_out['lonCell'] = IM['lonCell']
        ds_out['lonEdge'] = IM['lonEdge']
        ds_out['lonVertex'] = IM['lonVertex']
        ds_out['meshDensity'] = IM['meshDensity']
        ds_out['nEdgesOnCell'] = IM['nEdgesOnCell']
        ds_out['nEdgesOnEdge'] = IM['nEdgesOnEdge']
        ds_out['verticesOnCell'] = IM['verticesOnCell']
        ds_out['verticesOnEdge'] = IM['verticesOnEdge']
        ds_out['weightsOnEdge'] = IM['weightsOnEdge']
        ds_out['xCell'] = IM['xCell']
        ds_out['xEdge'] = IM['xEdge']
        ds_out['xVertex'] = IM['xVertex']
        ds_out['yCell'] = IM['yCell']
        ds_out['yEdge'] = IM['yEdge']
        ds_out['yVertex'] = IM['yVertex']
        ds_out['zCell'] = IM['zCell']
        ds_out['zEdge'] = IM['zEdge']
        ds_out['zVertex'] = IM['zVertex']
        
        ds_out['validOceanMask'] = mask
        ds_out['ismip6shelfMelt_3dThermalForcing'] = interptf
        ds_out['ismip6shelfMelt_zOcean'] = zlayers
        if (self.newLandIceFWFlux != 0):
            ds_out['floatingBasalMassBal'] = interpDS['floatingBasalMassBal'][:,:]
       
        # Save xtime
        #ds_out = ds_out.expand_dims({'StrLen': self.newXtime}, axis=1)

        print("self.newXtime: {}".format(self.newXtime))
        #ds_out = ds_out.expand_dims({'StrLen': self.newXtime})
        xtime = xr.DataArray(np.array(self.newXtime, dtype = np.dtype('S64')), dims=["Time"])
        xtime = xtime.encoding.update({"char_dim_name":"StrLen"})
        ds_out['xtime'] = xtime

        #clean history for new file
        if 'history' in ds_out.attrs:
            del ds_out.attrs['history']

        ds_out.to_netcdf(self.options.outputFile, mode='w', unlimited_dims=['Time'])
        ds_out.close()

        subprocess.run(["ncatted", "-a", "_FillValue,,d,,", self.options.outputFile])

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


