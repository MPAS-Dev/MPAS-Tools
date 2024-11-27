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
from argparse import ArgumentParser
import time
import cftime

class mpasToMaliInterp:
    
    def __init__(self):
        print("Gathering Information ... ")
        parser = ArgumentParser(
                prog='interpolate_mpasOcean_to_mali.py',
                description='Interpolates between MPAS-Ocean data to melt rates and thermal forcing on a MALI mesh')
        parser.add_argument("--oceanMesh", dest="mpasoMeshFile", help="MPAS-Ocean mesh to be interpolated from", metavar="FILENAME")
        parser.add_argument("--oceanDiags", dest="mpasoHistDir", help="Directory where MPAS-Ocean history files are stored. Should be only netcdf files in that directory", metavar="FILENAME")
        parser.add_argument("--ice", dest="maliFile", help="MALI mesh to be interpolated onto", metavar="FILENAME")
        parser.add_argument("-n", "--ntasks", dest="ntasks", help="Number of processors to use with ESMF_RegridWeightGen")
        parser.add_argument("-m", "--method", dest="method", help="Remapping method, either 'bilinear' or 'neareststod'")
        parser.add_argument("-o", "--outFile", dest="outputFile", help="Desired name of output file", metavar="FILENAME", default="mpas_to_mali_remapped.nc")
        parser.add_argument("-a","--yearlyAvg", dest="yearlyAvg", help="true or false option to average monthly output to yearly intervals", default="true")
        parser.add_argument("-s","--startYr", dest="startYr", type=int, help="starting year to process")
        parser.add_argument("-e","--endYr", dest="endYr", type=int, help="ending year to process (inclusive)")
        parser.add_argument("--meshVars", dest="includeMeshVars", help="whether to include mesh variables in resulting MALI forcing file", action='store_true')
        args = parser.parse_args()
        self.options = args

        self.mapping_file_name = f'mpaso_to_mali_mapping_{self.options.method}.nc'

        # Define ocean vertical coordinates for mali grid
        # For now, hardcode in depths from ISMIP6 AIS. Should work fine for our purposes
        self.ismip6shelfMelt_zOcean = np.array((-30, -90, -150, -210, -270, -330, -390, -450, -510,\
        -570, -630, -690, -750, -810, -870, -930, -990, -1050, -1110, -1170,\
        -1230, -1290, -1350, -1410, -1470, -1530, -1590, -1650, -1710, -1770), dtype='float64')

        # create bnds array
        self.ismip6shelfMelt_zBndsOcean = np.zeros((len(self.ismip6shelfMelt_zOcean), 2))
        self.ismip6shelfMelt_zBndsOcean[0, 0] = 0.0
        for z in range(1, len(self.ismip6shelfMelt_zOcean)):
            self.ismip6shelfMelt_zBndsOcean[z, 0] = 0.5 * (self.ismip6shelfMelt_zOcean[z - 1] +
                                                           self.ismip6shelfMelt_zOcean[z])
        for z in range(0, len(self.ismip6shelfMelt_zOcean) - 1):
            self.ismip6shelfMelt_zBndsOcean[z, 1] = 0.5 * (self.ismip6shelfMelt_zOcean[z] +
                                                           self.ismip6shelfMelt_zOcean[z + 1])
        self.ismip6shelfMelt_zBndsOcean[-1, 1] = self.ismip6shelfMelt_zOcean[-1] + \
                (self.ismip6shelfMelt_zOcean[-1] - self.ismip6shelfMelt_zBndsOcean[-1, 0])

    def prepare_mpaso_mesh_data(self):

        #open mpas ocean mesh
        OM = xr.open_dataset(self.options.mpasoMeshFile, decode_times=False, decode_cf=False)
        self.stTime = OM['simulationStartTime'].data.tobytes().decode()
        print(f'Using simulation start time of: {self.stTime}')
        self.landIceFloatingMask = OM['landIceFloatingMask'][0,:].data #not letting floating ice mask evolve for now because it's only in the mesh file

        # variables needed from MPAS-Ocean mesh file
        self.bottomDepth = OM['bottomDepth']
        self.maxLevelCell = OM['maxLevelCell']
        if 'minLevelCell' in OM:
           self.minLevelCell = OM['minLevelCell'][:].values
        else:
           self.minLevelCell = self.maxLevelCell * 0 + 1
        self.landIcePressure = OM['landIcePressure'][0,:].values

        self.nCells = OM.sizes['nCells']

        # create needed masks and remap them in separate file
        mpasoDomainMask = np.ones((self.nCells,), dtype=np.double)
        mpasoDomainMaskDA = xr.DataArray(mpasoDomainMask, name='mpasoDomainMask', dims=("nCells"), attrs={'long_name':'MPAS-Ocean domain mask'})
        mask1d = np.logical_not(self.landIceFloatingMask).astype(np.double)
        mpasoOpenOceanMaskDA = xr.DataArray(mask1d,
                                            name='mpasoOpenOceanMask',
                                            dims=("nCells",),
                                            attrs={'long_name':'MPAS-Ocean open ocean mask'})

        # create mask of 3d ocean data (including cavities)
        # Note: Technically, the 3d ocean mask will vary with time because fluctuations in
        # the layerThickess and surface pressure will cause the depths that are valid to
        # change slightly over time.  Eventually, we may want orig3dOceanMask to vary every
        # time step, but for now it is assumed to be constant in time.  To support the mask being
        # constant in time while the depths of the TF data is not, we have to ensure there is
        # valid TF data if any mismatch occurs (i.e. where the mask says valid data exists but
        # due to changes in layerThickness, it is outside the depth range of valid TF data).
        # This is handled by changing the vertical interpolation function to extrapolate instead
        # of insert nan when TF is calculated.  But here where we calculate the mask we mark
        # regions outside the valid depth range with nan so they can be masked.
        layerThickness = OM['layerThickness']
        mpas_cellCenterElev = compute_zmid(self.bottomDepth, self.maxLevelCell, layerThickness)
        vertInterpResult = _vertical_interpolate(self, np.zeros(layerThickness.shape), mpas_cellCenterElev.data, markExtrap=True)
        # the above commands follow the operations used for actual TF data below
        # the result will have zero where there is valid data and nan where there is not
        maskOcean3d = np.logical_not(np.isnan(vertInterpResult)).astype(np.double)
        maskTmp = np.swapaxes(maskOcean3d[0,:,:], 0, 1)  # necessary for ncremap
        mpaso3dOceanMaskDA = xr.DataArray(maskTmp,
                                          name='mpaso3dOceanMask',
                                          #dims=("nCells", "nISMIP6OceanLayers"),
                                          dims=("nISMIP6OceanLayers", "nCells"),
                                          attrs={'long_name':'MPAS-Ocean 3d ocean mask'})

        # prepare file to be remapped
        out_data_vars = xr.merge([mpasoDomainMaskDA, mpasoOpenOceanMaskDA, mpaso3dOceanMaskDA])
        dataOut = xr.Dataset(data_vars=out_data_vars)
        mpaso_mask_file = 'mpaso_mask_file.nc'
        dataOut.to_netcdf(mpaso_mask_file, mode='w')

        print("Calling ncremap for masks...")
        # remap the input data
        args_remap = ["ncremap",
                "-i", mpaso_mask_file,
                "-o", "mpaso_masks_on_mali_mesh.nc",
                "-m", self.mapping_file_name]
        check_call(args_remap)

        # Now on MALI mesh, clean up open ocean masks so they can be used by MALI extrap code
        ds = xr.open_dataset("mpaso_masks_on_mali_mesh.nc", decode_times=False, decode_cf=False)
        mask1d = ds.mpasoOpenOceanMask[:].values
        mask1d = np.where(mask1d > 0.99, 1.0, 0.0).astype(np.int32)  # only keep values close to 1
        mask2d = np.tile(mask1d.reshape(-1, 1), (1, len(self.ismip6shelfMelt_zOcean)))  # tile mask across all vertical layers
        mask3d = mask2d[np.newaxis, :, :]  # add time dimension
        mpasoOpenOceanMaskDA = xr.DataArray(mask3d,
                                            name='orig3dOceanMask',
                                            dims=("Time", "nCells", "nISMIP6OceanLayers"),
                                            attrs={'long_name':'MPAS-Ocean open ocean mask'})
        mpasoOpenOceanMaskDA.to_netcdf('orig3dOceanMask_open_ocean_mask.nc', mode='w')

        mask = ds.mpaso3dOceanMask[:].values
        mask = np.swapaxes(mask, 0, 1)  # undo swap prior to remapping
        mask = np.where(mask > 0.99, 1.0, 0.0).astype(np.int32)  # only keep values close to 1
        mask3d = mask[np.newaxis, :, :]  # add time dimension
        mpaso3dOceanMaskDA = xr.DataArray(mask3d,
                                          name='orig3dOceanMask',
                                          dims=("Time", "nCells", "nISMIP6OceanLayers"),
                                          attrs={'long_name':'MPAS-Ocean full 3d ocean mask'})
        mpaso3dOceanMaskDA.to_netcdf('orig3dOceanMask_3d_ocean_mask_with_cavities.nc', mode='w')


    def get_data(self, year):

        # open and concatenate diagnostic dataset
        self.DS = xr.open_mfdataset(os.path.join(self.options.mpasoHistDir,
            f'*.mpaso.hist.am.timeSeriesStatsMonthly.{year}-*-01.nc'), combine='nested', concat_dim='Time', decode_timedelta=False)
        
        # variables for time averaging
        self.temperature = self.DS['timeMonthly_avg_activeTracers_temperature']
        self.salinity = self.DS['timeMonthly_avg_activeTracers_salinity']
        self.density = self.DS['timeMonthly_avg_density']
        self.atmPressure =self. DS['timeMonthly_avg_atmosphericPressure']
        self.daysSinceStart = self.DS['timeMonthly_avg_daysSinceStartOfSim']
        avg_layerThickness = self.DS['timeMonthly_avg_layerThickness']
        self.layerThickness = avg_layerThickness.data

        #xt = DS['xtime_startMonthly']
        
        #xtime = np.array([xt],dtype=('S',64)) 
        # << NOTE >>: may need to eventually use: "xtime.data.tobytes().decode()" but not sure yet
        
        self.have_landIceFreshwaterFlux = False
        try:
            self.landIceFreshwaterFlux = self.DS['timeMonthly_avg_landIceFreshwaterFlux']
            self.have_landIceFreshwaterFlux = True
        except KeyError:
            print("No landIceFreshwaterFlux variable")

        self.coeff_0_cavity = self.DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_0']
        self.coeff_S_cavity = self.DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_S']
        self.coeff_p_cavity = self.DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_p']
        self.coeff_pS_cavity = self.DS.attrs['config_land_ice_cavity_freezing_temperature_coeff_pS']
        self.coeff_mushy_cavity = 0.0

        # Define vertical coordinates of mpas output
        mpas_cellCenterElev = compute_zmid(self.bottomDepth, self.maxLevelCell, avg_layerThickness)
        self.mpas_cellCenterElev = mpas_cellCenterElev.data

    def time_average_output(self):
        print("Time Averaging ...")

        #prepare time vector
        yearsSinceStart = self.daysSinceStart.values / 365.0
        stTime = np.datetime64(self.stTime[0:4])
        stYear = np.datetime_as_string(stTime,unit='Y').astype('float64')

        if (self.options.yearlyAvg == 'true'):

            #pre-allocate
            nt = 1
            _,nc,nz = self.temperature.shape
            self.newTemp = np.zeros((nt,nc,nz)) 
            self.newSal = np.zeros((nt,nc,nz))
            self.newDens = np.zeros((nt,nc,nz))
            self.newLThick= np.zeros((nt,nc,nz))
            self.newAtmPr = np.zeros((nt,nc))
            self.newMpasCCE = np.zeros((nt,nc,nz))
            if self.have_landIceFreshwaterFlux:
                self.newLandIceFWFlux = np.zeros((nt,nc))
            yearVec = np.zeros((nt,))    
            
            st = time.time()

            self.newTemp[0,:,:] = np.mean(self.temperature.values, axis=0)
            self.newSal[0,:,:] = np.mean(self.salinity.values, axis=0)
            self.newDens[0,:,:] = np.mean(self.density.values, axis=0)
            self.newLThick[0,:,:] = np.mean(self.layerThickness, axis=0)
            self.newMpasCCE[0,:,:] = np.mean(self.mpas_cellCenterElev, axis=0)
            self.newAtmPr[0,:] = np.mean(self.atmPressure.values, axis=0)
            if self.have_landIceFreshwaterFlux:
                self.newLandIceFWFlux[0,:] = np.mean(self.landIceFreshwaterFlux.values, axis=0)
            
            #Define time at the first of each year
            yearVec[0] = np.floor(np.min(yearsSinceStart)) + stYear

            nd = time.time()
            tm = nd - st
            print("Time averaging loop", tm, "seconds")
            
            #establish xtime
            dates = []
            xt = cftime.num2date(yearVec*365, units="days since 0000-01-01_00:00:00",calendar='noleap')
            xt_reformat = [dt.strftime("%Y-%m-%d_%H:%M:%S").ljust(64) for dt in xt]
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
            if self.have_landIceFreshwaterFlux:
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
                    kmax = self.maxLevelCell[iCell].values - 1

                    self.newPr[iTime,iCell,kmin] = self.newAtmPr[iTime,iCell] + \
                            self.landIcePressure[iCell] + \
                            self.newDens[iTime,iCell,kmin]*gravity*0.5*self.newLThick[iTime,iCell,kmin]

                    for k in np.arange(kmin + 1, kmax):
                        self.newPr[iTime,iCell,k] = self.newPr[iTime,iCell,k-1] + \
                                0.5 * gravity * (self.newDens[iTime,iCell,k-1] * self.newLThick[iTime,iCell,k-1] \
                                + self.newDens[iTime,iCell,k] * self.newLThick[iTime,iCell,k])

                    ocn_freezing_temperature[iTime,iCell,:] = self.coeff_0_cavity + \
                            self.coeff_S_cavity * self.newSal[iTime,iCell,:] + \
                            self.coeff_p_cavity * self.newPr[iTime,iCell,:] \
                            + self.coeff_pS_cavity * self.newPr[iTime,iCell,:] * self.newSal[iTime,iCell,:] \
                            + self.coeff_mushy_cavity * self.newSal[iTime,iCell,:] / (1.0 - self.newSal[iTime,iCell,:] / 1e3)
        nd = time.time()
        tm = nd - st
        print("Ocean thermal forcing loop", tm, "seconds")
        # Calculate thermal forcing
        self.oceanThermalForcing = self.newTemp - ocn_freezing_temperature

    def create_mapping_file(self):

        #create scrip files
        print("Creating Scrip Files")
        mpaso_scripfile = "tmp_mpaso_scrip.nc"
        mali_scripfile = "tmp_mali_scrip.nc"

        scrip_from_mpas(self.options.mpasoMeshFile, mpaso_scripfile)
        scrip_from_mpas(self.options.maliFile, mali_scripfile)

        #create mapping file
        print("Creating Mapping File")
        args_esmf = ['srun',
                '-n', self.options.ntasks, 'ESMF_RegridWeightGen',
                '--source', mpaso_scripfile,
                '--destination', mali_scripfile,
                '--weight', self.mapping_file_name,
                '--method', self.options.method,
                '-i', '-64bit_offset',
                "--dst_regional", "--src_regional", '--ignore_unmapped']
        check_call(args_esmf)
    
    def remap_mpas_to_mali(self, year):

        out_name = f'{self.options.outputFile}_{year}.nc'
        tmp_mpasoSourceFile = "tmp_mpasoSourceFile.nc"
        tmp_maliDestFile = "tmp_maliDestFile.nc"

        print("Start remapping ... ")

        # Get data needed for remapping
        f = xr.open_dataset(self.options.mpasoMeshFile, decode_times=False, decode_cf=False)
        tf = xr.DataArray(self.oceanThermalForcing.astype('float64'),dims=('Time','nCells','nVertLevels'))
        mcce = xr.DataArray(self.newMpasCCE.astype('float64'),dims=('Time','nCells','nVertLevels'))

        # perform vertical interpolation
        # do this before horiz interp to avoid any horiz/vert "mixing" due to
        # potential mpas-ocean hybrid coordinate.  Also should speed up ncremap
        # because the number of ismip6 vert levels is less than most ocean meshes
        st = time.time()
        vertInterpTF = _vertical_interpolate(self, tf.data, mcce.data)
        nd = time.time()
        tm = nd - st
        print("vertical interpolation:", tm, "seconds")

        st = time.time()
        # create new file of ocean data to be remapped
        tf_DA = xr.DataArray(vertInterpTF.astype('float64'),
                             dims=('Time','nCells','nISMIP6OceanLayers'))
        dataOut = xr.Dataset(data_vars={'ismip6shelfMelt_3dThermalForcing': tf_DA})
        if self.have_landIceFreshwaterFlux:
            melt = xr.DataArray(self.newLandIceFWFlux.astype('float64'),dims=('Time','nCells'))
            dataOut['floatingBasalMassBal'] = -1.0 * melt
        dataOut.to_netcdf(tmp_mpasoSourceFile, mode='w')
        f.close()

        subprocess.run(["ncatted", "-a", "_FillValue,,d,,", tmp_mpasoSourceFile])
        subprocess.run(["ncpdq", "-O", "-a", "Time,nISMIP6OceanLayers,nCells", tmp_mpasoSourceFile, tmp_mpasoSourceFile])

        nd = time.time()
        tm = nd - st
        print(f"saved modified mpaso data file for remapping and file prep: {tm} seconds")

        print("ncremap ...")
        st = time.time()
        # remap the input data
        args_remap = ["ncremap",
                "-i", tmp_mpasoSourceFile,
                "-o", tmp_maliDestFile,
                "-m", self.mapping_file_name]
        check_call(args_remap)
        subprocess.run(["ncpdq", "-O", "-a", "Time,nCells,nISMIP6OceanLayers", tmp_maliDestFile, tmp_maliDestFile])
        nd = time.time()
        tm = nd - st
        print(f"ncremap completed: {tm} seconds")

        print(f"Saving data to {out_name}")

        #create new mali forcing file that only contains desired forcing variables
        ds_remapped = xr.open_dataset(tmp_maliDestFile, decode_times=False, decode_cf=False)
        ds_out = xr.Dataset(data_vars={'ismip6shelfMelt_3dThermalForcing': ds_remapped['ismip6shelfMelt_3dThermalForcing']})
        ds_out['ismip6shelfMelt_zOcean'] = xr.DataArray(self.ismip6shelfMelt_zOcean, dims=("nISMIP6OceanLayers"))
        ds_out['ismip6shelfMelt_zBndsOcean'] = xr.DataArray(self.ismip6shelfMelt_zBndsOcean, dims=("nISMIP6OceanLayers", "TWO"))
        if self.have_landIceFreshwaterFlux:
            ds_out['floatingBasalMassBal'] = ds_remapped['floatingBasalMassBal']

        if self.options.includeMeshVars:
            IM = xr.open_dataset(self.options.maliFile,decode_times=False,decode_cf=False)
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
            IM.close()
        
        # Save xtime
        #ds_out = ds_out.expand_dims({'StrLen': self.newXtime}, axis=1)

        print("self.newXtime: {}".format(self.newXtime))
        #ds_out = ds_out.expand_dims({'StrLen': self.newXtime})
        xtime = xr.DataArray(np.array(self.newXtime, dtype = np.dtype('S64')), dims=["Time"])
        xtime.encoding.update({"char_dim_name":"StrLen"})
        ds_out['xtime'] = xtime

        ds_out.to_netcdf(out_name, mode='w', unlimited_dims=['Time'])
        ds_out.close()

        subprocess.run(["ncatted", "-a", "_FillValue,,d,,", out_name])

        #remove temporary files
        files = os.listdir()
        for file in files:
            if file.startswith("tmp"):
                os.remove(file)

def _vertical_interpolate(self, TF, mpas_cellCenterElev, markExtrap=False):
        #Vertical interpolation
        print("Vertically interpolating onto standardized z-level grid")

        nt,nc,nz1 = mpas_cellCenterElev.shape
        nz2, = self.ismip6shelfMelt_zOcean.shape

        # numpy.interp expects the x-coordinate sequence to be increasing
        # This is not explicitly enforced, and it's unclear if it is
        # necessary, but we will follow that convention.
        # So we need to flip the ordering of the zOcean array
        z_target = np.flip(self.ismip6shelfMelt_zOcean, 0)
        assert(np.all(np.diff(z_target) > 0))

        # Reshape to (nt*nc, nz) before interpolation to avoid looping
        # Also, flip the ordering of the z-dimension
        # numpy.interp requires the independent variable to be increasing
        # but cellCenterElev will be decreasing (it is indexed from ocean
        # surface down with positive up)
        TF = TF.reshape(nt*nc, -nz1)
        TF = np.flip(TF, 1)
        mpas_cellCenterElev = mpas_cellCenterElev.reshape(nt*nc, -1)
        mpas_cellCenterElev = np.flip(mpas_cellCenterElev, 1)

        # Linear interpolation
        nan_count = 0
        vertInterpTF = np.zeros((nt*nc, nz2))
        for i in range(nt*nc):
            ind = np.where(~np.isnan(TF[i,:].flatten() * mpas_cellCenterElev[i,:].flatten()))

            if len(ind[0]) != 0:
                assert(np.all(np.diff(mpas_cellCenterElev[i,ind]) > 0))  # confirm correct ordering
                # Note: Setting locations outside the extent of valid ocean data to nan
                # so that those locations can be avoided
                # This linear interpolation approach means that any vertical levels on the destination
                # vertical grid that are not bounded on both sides by valid vertical levels on the
                # source grid will be marked as nan.  We may want to replace this method with one
                # that uses and source data overlapping the range of the destination vertical level,
                # as well as doing an area weighted remapping of the overlap.
                if markExtrap:
                    vertInterpTF[i,:] = np.interp(z_target, mpas_cellCenterElev[i,ind].flatten(), TF[i,ind].flatten(), right=np.nan, left=np.nan)
                else:
                    vertInterpTF[i,:] = np.interp(z_target, mpas_cellCenterElev[i,ind].flatten(), TF[i,ind].flatten())
            else:
                nan_count = nan_count + 1

        # Reshape vertInterpTF back to (nt, nc, nz), flip back vertical coordinate
        vertInterpTF = vertInterpTF.reshape(nt, nc, nz2)
        vertInterpTF = np.flip(vertInterpTF, 2)

        return vertInterpTF

def main():
        run = mpasToMaliInterp()

        run.create_mapping_file()

        run.prepare_mpaso_mesh_data()

        for yr in range(run.options.startYr, run.options.endYr+1):
           st = time.time()
           print(f'\n**** Processing year {yr} ****\n')

           run.get_data(yr)
        
           #compute yearly time average if necessary
           run.time_average_output()
  
           #calculate thermal forcing
           run.calc_ocean_thermal_forcing()

           #remap mpas to mali
           run.remap_mpas_to_mali(yr)

           run.DS.close()


           nd = time.time()
           tm = nd-st
           print(f"\n### Processed year {yr} in {tm} seconds\n")

        print("Finished all years.")
if __name__ == "__main__":
    main()


