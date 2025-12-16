#!/usr/bin/env python

import os
import shutil
import xarray as xr
import subprocess
import numpy as np
import pandas as pd
from scipy import stats as st
from mpas_tools.logging import check_call
from mpas_tools.scrip.from_mpas import scrip_from_mpas
from mpas_tools.ocean.depth import compute_zmid
from argparse import ArgumentParser
import time
import cftime
import glob
import json
from shapely.geometry import shape
from shapely.prepared import prep
from shapely.geometry import Point

class eccoToMaliInterp:
    
    def __init__(self):
        print("Gathering Information ... ")
        parser = ArgumentParser(
                prog='interpolate_ecco_to_mali.py',
                description='Interpolates ECCO reanalysis data onto a MALI mesh')
        parser.add_argument("--eccoDir", dest="eccoDir", required=True, help="Directory where ECCO files are stored. Should be only netcdf files in that directory", metavar="FILENAME")
        parser.add_argument("--maliMesh", dest="maliMesh", required=True, help="MALI mesh to be interpolated onto", metavar="FILENAME")
        parser.add_argument("--eccoVars", type=str, nargs='+', dest="eccoVars", required=True, help="ECCO variables to interpolate, current options are any combination of 'THETA', 'SALT', 'SIarea'")
        parser.add_argument("--geojson", dest="geojson", help="Optional geojson used to permanently define icebergFjordMask within a region")
        parser.add_argument("--maxDepth", type=int, dest="maxDepth", default=1500, help="Maximum depth in meters to include in the interpolation")
        parser.add_argument("--ntasks", dest="ntasks", type=str, help="Number of processors to use with ESMF_RegridWeightGen", default='128')
        parser.add_argument("--method", dest="method", type=str, help="Remapping method, either 'bilinear' or 'neareststod'", default='conserve')
        parser.add_argument("--outFile", dest="outputFile", help="Desired name of output file", metavar="FILENAME", default="ecco_to_mali_remapped.nc")
        parser.add_argument("--meshVars", dest="includeMeshVars", help="whether to include mesh variables in resulting MALI forcing file", action='store_true')
        parser.add_argument("--iceAreaThresh", dest="iceAreaThresh", type=float, help="Threshold sea ice fraction used to define icebergFjordMask.", default=0.5)
        parser.add_argument("--extrapIters", dest="extrapIters", type=str, help="Number of sea ice extrapolation iterations when regridding", default='50')
        parser.add_argument("--keepTmpFiles", dest="keepTmpFiles", help="Whether to keep temporary netcdf files created throughout script. Useful for debugging", action="store_true")
        args = parser.parse_args()
        self.options = args

        self.mapping_file_name = f'ecco_to_mali_mapping_{self.options.method}.nc'
        self.tileNums = [14, 27] # Hardcoding ECCO tiles for now. Indexing is dependent on tile order, do not change

    def merge_MITgcm_files_to_unstruct_grid(self):
        ncfiles = glob.glob(os.path.join(self.options.eccoDir, "*.nc"))
        saltBool = 0
        thetaBool = 0
        siaBool = 0
    
        #loop through MITgcm tile files and create unstructured arrays for each variable in eccoVars
        for tile in self.tileNums:
            print(f"Processing Tile {tile}")
            o3dmBool = 0
            gridBool = 0
            for file in ncfiles:
                print(f"Processing File {file}")
                paddedTile = str(tile).zfill(4)
                print(f"Padded Tile: {paddedTile}")
                if paddedTile in file:
                    print("Padded tile is in file")
                    
                    if 'SALT' in file and 'SALT' in self.options.eccoVars:
                        validFile = self.options.eccoDir + '/SALT.' + paddedTile + '.nc'
                        ds_ecco = xr.open_dataset(validFile)
                        
                        # only open depth variable once, regardless of which variables we are interested in 
                        if saltBool == 0 and thetaBool == 0:
                            z = ds_ecco['dep'].values
                            
                        if tile == 14:
                            sal = ds_ecco['SALT'][:].values
                            # omit rows/columns where no corner information in scrip file
                            sal = sal[:,:,:,1:-1]
                            nt,nz,nx,ny = np.shape(sal)
                            sal = sal.reshape(nt, nz, nx * ny)
                        elif tile == 27:
                            sl = ds_ecco['SALT'][:].values
                            # omit rows/columns where no corner information in scrip file
                            sl = sl[:,:,0:-1,0:-1]
                            nt,nz,nx,ny = np.shape(sl)
                            sl = sl.reshape(nt, nz, nx * ny)
                            sal = np.concatenate((sal,sl), axis=2)
                        
                        # ECCO 'land' variable is equivalent to MALI's orig3dOceanMask, despite the counterintuitive name.
                        # Only need to call this once per tile, but needs to be using a theta/salinity file because it includes depth 
                        # No need to call this if only interpolating SIarea
                        if tile == 14 and o3dmBool == 0:
                            orig3dOceanMask = ds_ecco['land'].values
                            # omit rows/columns where no corner information in scrip file
                            orig3dOceanMask = orig3dOceanMask[:,:,1:-1]
                            nz,nx,ny = np.shape(orig3dOceanMask)
                            orig3dOceanMask = orig3dOceanMask.reshape(nz, nx * ny)
                            od3mBool = 1
                        elif tile == 27 and o3dmBool == 0:
                            o3dm = ds_ecco['land'].values
                            # omit rows/columns where no corner information in scrip file
                            o3dm = o3dm[:,0:-1,0:-1]
                            nz,nx,ny = np.shape(o3dm)
                            o3dm = o3dm.reshape(nz, nx * ny)
                            orig3dOceanMask = np.concatenate((orig3dOceanMask, o3dm), axis=1)
                            o3dmBool = 1
                        
                        # set MALI invalid value
                        boolZ,boolXY = np.where(orig3dOceanMask == 0)
                        for iZ,iXY in zip(boolZ,boolXY):
                            sal[:,iZ,iXY] = 1e36
                        
                        saltBool = 1       
                    
                    # Repeat for potential temperature if needed
                    elif 'THETA' in file and 'THETA' in self.options.eccoVars:       
                        validFile = self.options.eccoDir + '/THETA.' + paddedTile + '.nc'
                        ds_ecco = xr.open_dataset(validFile)
                        
                        if saltBool == 0 and thetaBool == 0:
                            z = ds_ecco['dep'].values
                            
                        if tile == 14:
                            theta = ds_ecco['THETA'][:].values
                            # omit rows/columns where no corner information in scrip file
                            theta = theta[:,:,:,1:-1]
                            nt,nz,nx,ny = np.shape(theta)
                            theta = theta.reshape(nt, nz, nx * ny)
                        elif tile == 27:
                            th = ds_ecco['THETA'][:].values
                            # omit rows/columns where no corner information in scrip file
                            th = th[:,:,0:-1,0:-1]
                            nt,nz,nx,ny = np.shape(th)
                            th = th.reshape(nt, nz, nx * ny)
                            theta = np.concatenate((theta,th), axis=2)

                        # ECCO 'land' variable is equivalent to MALI's orig3dOceanMask, despite the counterintuitive name.
                        # Only need to call this once per tile, but needs to be using a theta/salinity file because it includes depth 
                        # No need to call this if only interpolating SIarea
                        if tile == 14 and o3dmBool == 0:
                            orig3dOceanMask = ds_ecco['land'].values
                            # omit rows/columns where no corner information in scrip file
                            orig3dOceanMask = orig3dOceanMask[:,:,1:-1]
                            nz,nx,ny = np.shape(orig3dOceanMask)
                            orig3dOceanMask = orig3dOceanMask.reshape(nz, nx * ny)
                            od3mBool = 1
                        elif tile == 27 and o3dmBool == 0:
                            o3dm = ds_ecco['land'].values
                            # omit rows/columns where no corner information in scrip file
                            o3dm = o3dm[:,0:-1,0:-1]
                            nz,nx,ny = np.shape(o3dm)
                            o3dm = o3dm.reshape(nz, nx * ny)
                            orig3dOceanMask = np.concatenate((orig3dOceanMask, o3dm), axis=1)
                            o3dmBool = 1
                        
                        # set MALI invalid value
                        boolZ,boolXY = np.where(orig3dOceanMask == 0)
                        for iZ,iXY in zip(boolZ,boolXY):
                            theta[:,iZ,iXY] = 1e36
                        
                        thetaBool = 1
                        
                    # Repeat for sea ice fraction if needed
                    elif 'SIarea' in file and 'SIarea' in self.options.eccoVars:
                        siaBool = 1
                        validFile = self.options.eccoDir + '/SIarea.' + paddedTile + '.nc'
                        ds_ecco = xr.open_dataset(validFile)
                        if tile == 14:
                            sia = ds_ecco['SIarea'][:].values
                            # omit rows/columns where no corner information in scrip file
                            sia = sia[:,:,1:-1]
                            nt,nx,ny = np.shape(sia)
                            sia = sia.reshape(nt, nx * ny)
                        elif tile == 27:
                            si = ds_ecco['SIarea'][:].values
                            # omit rows/columns where no corner information in scrip file
                            si = si[:,0:-1,0:-1]
                            nt,nx,ny = np.shape(si)
                            si = si.reshape(nt, nx * ny)
                            sia = np.concatenate((sia,si), axis=1)

            if ('SALT' in self.options.eccoVars and saltBool == 0):
                raise ValueError(f"SALT netcdf file for tile {tile} not found in {self.options.eccoDir}")

            if ('THETA' in self.options.eccoVars and thetaBool == 0):
                raise ValueError(f"THETA netcdf file for tile {tile} not found in {self.options.eccoDir}")
                
            if ('SIarea' in self.options.eccoVars and siaBool == 0):
                raise ValueError(f"SIarea netcdf file for tile {tile} not found in {self.options.eccoDir}")

        # combine into new netcdf. Need to load depth and time from original ecco files. Assuming time/depth is consistent
        # across all files, just load info from most recent
        
        # Establish dataset, and fill in standard dimensions. other dimensions are added on an as-needed basis
        ds_unstruct = xr.Dataset()
        ds_unstruct.expand_dims(["nCells", "Time"])
        
        # Save ECCO time variable as datetime string. Will be converted into xtime format when saving final output file
        time = ds_ecco['tim'].values
        self.xtime_str = [pd.to_datetime(dt).strftime("%Y-%m-%d_%H:%M:%S").ljust(64) for dt in time]
        
        # save unstructured variables with MALI/MPAS-O names 
        if 'z' in locals():
            ds_unstruct.expand_dims("nISMIP6OceanLayers")
            ds_unstruct.expand_dims("TWO")
            z = z[z <= self.options.maxDepth]
        
            # create bnds array
            zBnds = np.zeros((len(z), 2))
            zBnds[0, 0] = 0.0
            for i in range(1, len(z)):
                zBnds[i, 0] = 0.5 * (z[i - 1] + z[i])
            for i in range(0, len(z) - 1):
                zBnds[i, 1] = 0.5 * (z[i] + z[i + 1])
            zBnds[-1, 1] = z[-1] + (z[-1] - zBnds[-1, 0])
            
            DA_z = xr.DataArray(-z.astype('float64'), dims=("nISMIP6OceanLayers")) #MALI wants negative depth coordinates
            DA_z.attrs['long_name'] = 'depth'
            DA_z.attrs['units'] = 'm'
            ds_unstruct['ismip6shelfMelt_zOcean'] = DA_z

            DA_zBnds = xr.DataArray(-zBnds.astype('float64'), dims=("nISMIP6OceanLayers", "TWO")) #MALI wants negative depth coordinates
            DA_zBnds.attrs['long_name'] = 'Bounds for ismip6 ocean layers'
            DA_zBnds.attrs['units'] = 'm'
            ds_unstruct['ismip6shelfMelt_zBndsOcean'] = DA_zBnds
        
        if 'theta' in locals():
            print(f"z shape: {z.shape}")
            print(f"theta shape: {theta.shape}")
            theta = theta[:,0:len(z),:]
            print(f"new theta shape: {theta.shape}")
            DA_theta = xr.DataArray(theta.astype('float64'), dims=("Time", "nISMIP6OceanLayers", "nCells"))
            DA_theta.attrs['long_name'] = 'Potential_Temperature'
            DA_theta.attrs['units'] = 'degC'
            ds_unstruct['oceanTemperature'] = DA_theta
                
        if 'sal' in locals():
            sal = sal[:,0:len(z),:]
            DA_salt = xr.DataArray(sal.astype('float64'), dims=("Time", "nISMIP6OceanLayers", "nCells"))
            DA_salt.attrs['long_name'] = 'Salinity'
            DA_salt.attrs['units'] = 'psu'
            ds_unstruct['oceanSalinity'] = DA_salt
                
        if 'sia' in locals():
            DA_sia = xr.DataArray(sia.astype('float64'), dims=("Time", "nCells"))
            DA_sia.attrs['long_name'] = 'Sea ice fractional ice-covered area [0 to 1]'
            DA_sia.attrs['units'] = 'm2/m2'
            # create separate netcdf file for SIA so it can be extrapolated during remapping
            ds_sia = xr.Dataset()
            ds_sia['iceAreaCell'] = DA_sia
            self.siaFile = 'sia.tmp.nc'
            ds_sia.to_netcdf(self.siaFile)
            ds_sia.close()

        if 'orig3dOceanMask' in locals():
            orig3dOceanMask = orig3dOceanMask[0:len(z),:]
            orig3dOceanMask = np.tile(orig3dOceanMask[np.newaxis, :, :], (nt, 1, 1))
            DA_o3dm = xr.DataArray(orig3dOceanMask.astype('float64'), dims=("Time", "nISMIP6OceanLayers", "nCells"))
            ds_unstruct['orig3dOceanMask'] = DA_o3dm 

        self.ecco_unstruct = "ecco_combined_unstructured.tmp.nc"
        ds_unstruct.to_netcdf(self.ecco_unstruct, unlimited_dims=['Time'])
        ds_unstruct.close()
        print("ECCO restructuring complete")
        
    def remap_ecco_to_mali(self):
        print("Creating MALI scrip file")
        mali_scripfile = "mali.scrip.nc"
        scrip_from_mpas(self.options.maliMesh, mali_scripfile)

        #create mapping file
        print("Creating Primary Mapping File")
        args = ['srun',
                '-n', self.options.ntasks, 'ESMF_RegridWeightGen',
                '--source', self.ecco_scripfile,
                '--destination', mali_scripfile,
                '--weight', self.mapping_file_name,
                '--method', self.options.method,
                '--netcdf4',
                "--dst_regional", "--src_regional", '--ignore_unmapped']
        check_call(args)
    
        #remap
        print("Start remapping ... ")

        subprocess.run(["ncatted", "-a", "_FillValue,,d,,", self.ecco_unstruct])
        subprocess.run(["ncpdq", "-O", "-a", "Time,nISMIP6OceanLayers,nCells", self.ecco_unstruct, self.ecco_unstruct])

        print("ncremap ...")
        st = time.time()
        # remap the input data
        remappedOutFile = self.options.outputFile[0:-3] + '.tmp.remapped.nc'
        args_remap = ["ncremap",
                "-i", self.ecco_unstruct,
                "-o", remappedOutFile,
                "-m", self.mapping_file_name]
        check_call(args_remap)

        subprocess.run(["ncpdq", "-O", "-a", "Time,nCells,nISMIP6OceanLayers", remappedOutFile, remappedOutFile])
        nd = time.time()
        tm = nd - st
        print(f"ncremap completed: {tm} seconds")
        print("Primary Remapping completed")

        # repeat remapping for sia, including extrapolation of sia to fill domain
        ds_siaScrip = xr.open_dataset(self.ecco_scripfile)
        ds_unstruct = xr.open_dataset(self.ecco_unstruct)
        o3dm = ds_unstruct['orig3dOceanMask'][1,1,:].values
        valid = (o3dm > 0.9999).astype('int32') # Account for rounding error
        ds_siaScrip['grid_imask'] = xr.DataArray(valid.astype('int32'),
                                                        dims=('grid_size'),
                                                        attrs={'units':'unitless'})
        sia_scripfile = 'ecco_sia.scrip.nc'
        ds_siaScrip.to_netcdf(sia_scripfile)
        ds_siaScrip.close()

        print("Creating SIA mapping File")
        sia_mapping_file = 'ecco_to_mali_mapping_bilinear_sia.nc'
        args = ['srun',
                '-n', self.options.ntasks, 'ESMF_RegridWeightGen',
                '--source', sia_scripfile,
                '--destination', mali_scripfile,
                '--weight', sia_mapping_file,
                '--method', 'bilinear',
                '--extrap_method', 'creep',
                '--extrap_num_levels', self.options.extrapIters,
                '--netcdf4',
                "--dst_regional", "--src_regional", '--ignore_unmapped']
        check_call(args)
    
        #remap
        print("Start SIA remapping ... ")

        subprocess.run(["ncatted", "-a", "_FillValue,,d,,", self.siaFile])

        print("ncremap of SIA...")
        st = time.time()
        # remap the input data
        siaRemappedOutFile = self.options.outputFile[0:-3] + '.tmp.sia.remapped.nc'
        args_remap = ["ncremap",
                "-i", self.siaFile,
                "-o", siaRemappedOutFile,
                "-m", sia_mapping_file]
        check_call(args_remap)

        nd = time.time()
        tm = nd - st
        print(f"ncremap completed: {tm} seconds")
        print("SIA Remapping completed")
        
        print("Applying final edits to remapped file ...")
        # define icebergFjordMask from SIA
        ds_out = xr.open_dataset(remappedOutFile, decode_times=False, decode_cf=False, mode='r+')
        ds_sia = xr.open_dataset(siaRemappedOutFile)
        iac = ds_sia['iceAreaCell'].values
        icebergFjordMask = np.zeros(iac.shape)
        # Find cells where sea ice fraction exceeds threshold
        seaice_mask = iac > self.options.iceAreaThresh 
        icebergFjordMask[seaice_mask] = 1

        # Use geojson to add any additional cells to icebergFjordMask
        if hasattr(self.options, "geojson"):
            with open(self.options.geojson) as f:
                geo = json.load(f)
                polygon = prep(shape(geo["features"][0]["geometry"]))
            
            lonCell = ds_out['lon'][:].values.ravel() - 360 # Convert to geojson longitude convention
            latCell = ds_out['lat'][:].values.ravel()
            
            points = [Point(xy) for xy in zip(lonCell, latCell)]
            geo_mask = np.array([polygon.covers(p) for p in points]).astype('bool')
            icebergFjordMask[:,geo_mask] = 1

        ifm = xr.DataArray(icebergFjordMask.astype('int32'), dims=("Time", "nCells"))
        sia = xr.DataArray(iac.astype('float32'), dims=("Time", "nCells"))

        # During conservative remapping, some valid ocean cells are averaged with invalid ocean cells. Use
        # float version of orig3dOceanMask to identify these cells and remove. Resave orig3dOceanMask as int32
        # so MALI recognizes it
        o3dm = ds_out['orig3dOceanMask']
        temp = ds_out['oceanTemperature']
        sal = ds_out['oceanSalinity']

        mask = o3dm > 0.9999 #Account for rounding error. Good cells are not exactly 1
        ds_out['orig3dOceanMask'] = xr.DataArray(mask.astype('int32'), dims=("Time", "nCells", "nISMIP6OceanLayers"))
        ds_out['oceanTemperature'] = temp.where(mask, 1e36).astype('float64')
        ds_out['oceanSalinity'] = sal.where(mask, 1e36).astype('float64')
        # Add sia and icebergBergMask into main file
        ds_out['iceAreaCell'] = sia
        ds_out['icebergFjordMask'] = ifm
        
        # save xtime correctly
        xtime = np.array(self.xtime_str, dtype = np.dtype(('S',64)))
        DA_xtime = xr.DataArray(xtime, dims=["Time"])
        DA_xtime.encoding.update({"char_dim_name": "StrLen"})
        ds_out['xtime'] = DA_xtime

        # Save final output file if not adding mesh variables
        if not self.options.includeMeshVars:
            # <<NOTE>>: Creating a new netcdf file like this changes the format of xtime. Need to add
            # encoding update to maintain original formatting
            ds_out.to_netcdf(self.options.outputFile, unlimited_dims=['Time'])
            subprocess.run(["ncatted", "-a", "_FillValue,,d,,", altRemappedOutFile])

        else:
            # Transfer MALI mesh variables if needed
            ds_mali = xr.open_dataset(self.options.maliMesh, decode_times=False, decode_cf=False)
            ds_out['angleEdge'] = ds_mali['angleEdge']
            ds_out['areaCell'] = ds_mali['areaCell']
            ds_out['areaTriangle'] = ds_mali['areaTriangle']
            ds_out['cellsOnCell'] = ds_mali['cellsOnCell']
            ds_out['cellsOnEdge'] = ds_mali['cellsOnEdge']
            ds_out['cellsOnVertex'] = ds_mali['cellsOnVertex']
            ds_out['dcEdge'] = ds_mali['dcEdge']
            ds_out['dvEdge'] = ds_mali['dvEdge']
            ds_out['edgesOnCell'] = ds_mali['edgesOnCell']
            ds_out['edgesOnEdge'] = ds_mali['edgesOnEdge']
            ds_out['edgesOnVertex'] = ds_mali['edgesOnVertex']
            ds_out['indexToCellID'] = ds_mali['indexToCellID']
            ds_out['indexToEdgeID'] = ds_mali['indexToEdgeID']
            ds_out['indexToVertexID'] = ds_mali['indexToVertexID']
            ds_out['kiteAreasOnVertex'] = ds_mali['kiteAreasOnVertex']
            ds_out['latCell'] = ds_mali['latCell']
            ds_out['latEdge'] = ds_mali['latEdge']
            ds_out['latVertex'] = ds_mali['latVertex']
            ds_out['lonCell'] = ds_mali['lonCell']
            ds_out['lonEdge'] = ds_mali['lonEdge']
            ds_out['lonVertex'] = ds_mali['lonVertex']
            ds_out['meshDensity'] = ds_mali['meshDensity']
            ds_out['nEdgesOnCell'] = ds_mali['nEdgesOnCell']
            ds_out['nEdgesOnEdge'] = ds_mali['nEdgesOnEdge']
            ds_out['verticesOnCell'] = ds_mali['verticesOnCell']
            ds_out['verticesOnEdge'] = ds_mali['verticesOnEdge']
            ds_out['weightsOnEdge'] = ds_mali['weightsOnEdge']
            ds_out['xCell'] = ds_mali['xCell']
            ds_out['xEdge'] = ds_mali['xEdge']
            ds_out['xVertex'] = ds_mali['xVertex']
            ds_out['yCell'] = ds_mali['yCell']
            ds_out['yEdge'] = ds_mali['yEdge']
            ds_out['yVertex'] = ds_mali['yVertex']
            ds_out['zCell'] = ds_mali['zCell']
            ds_out['zEdge'] = ds_mali['zEdge']
            ds_out['zVertex'] = ds_mali['zVertex']
            # <<NOTE>>: Creating a new netcdf file like this changes the format of xtime. Need to add
            # encoding update to maintain original formatting
            ds_out.to_netcdf(self.options.outputFile, unlimited_dims=['Time'])
            ds_out.close()
            ds_mali.close()
            subprocess.run(["ncatted", "-a", "_FillValue,,d,,", self.options.outputFile])
        
        # Remove temporary output files
        if not self.options.keepTmpFiles :
            for filename in os.listdir("."):
                if "tmp" in filename and os.path.isfile(filename):
                    os.remove(filename)

    def create_ECCO_scrip_file(self):
        
        print("Creating ECCO scrip file")

        for tile in self.tileNums :
            # start processing with grid files
            paddedTile = str(tile).zfill(4)
            gridFile = self.options.eccoDir + '/GRID.' + paddedTile + '.nc'
            grd = xr.open_dataset(gridFile)

            # identify cell centers and corners
            # lat/lon centers
            XC = grd['XC'][:].values
            YC = grd['YC'][:].values
            
            # lat/lon corners
            XG = grd['XG'][:].values
            YG = grd['YG'][:].values

            Xse_corner = np.zeros((XC.shape))
            Xsw_corner = np.zeros((XC.shape))
            Xnw_corner = np.zeros((XC.shape))
            Xne_corner = np.zeros((XC.shape))

            Yse_corner = np.zeros((YC.shape))
            Ysw_corner = np.zeros((YC.shape))
            Ynw_corner = np.zeros((YC.shape))
            Yne_corner = np.zeros((YC.shape))

            nx,ny = XC.shape
            for i in np.arange(0,nx-1):
                for ii in np.arange(0,ny-1):
                    Xse_corner[i,ii] = XG[i,ii]
                    Xsw_corner[i,ii] = XG[i+1,ii]
                    Xnw_corner[i,ii] = XG[i+1,ii+1]
                    Xne_corner[i,ii] = XG[i,ii+1]

                    Yse_corner[i,ii] = YG[i,ii]
                    Ysw_corner[i,ii] = YG[i+1,ii]
                    Ynw_corner[i,ii] = YG[i+1,ii+1]
                    Yne_corner[i,ii] = YG[i,ii+1]
                    
            if tile == 14:
                # Stitch together corner info from neighboring tile
                gridFile = self.options.eccoDir + '/GRID.0027.nc'
                grd27 = xr.open_dataset(gridFile)
                xg27 = grd27['XG'][:].values
                yg27 = grd27['YG'][:].values
                xg = np.flip(xg27[1:,0])
                yg = np.flip(yg27[1:,0])

                for i in np.arange(0,ny-1):        
                    Xse_corner[nx-1,i] = XG[nx-1,i]
                    Xne_corner[nx-1,i] = XG[nx-1,i+1]
                    Xnw_corner[nx-1,i] = xg[i]
                    Xsw_corner[nx-1,i] = xg[i-1]
                    
                    Yse_corner[nx-1,i] = YG[nx-1,i]
                    Yne_corner[nx-1,i] = YG[nx-1,i+1] 
                    Ynw_corner[nx-1,i] = yg[i]
                    Ysw_corner[nx-1,i] = yg[i-1]
        
                # Remove edge cells without corner information
                XC = XC[:,1:-1]
                YC = YC[:,1:-1]
                Xse_corner = Xse_corner[:,1:-1]
                Xne_corner = Xne_corner[:,1:-1]
                Xnw_corner = Xnw_corner[:,1:-1]
                Xsw_corner = Xsw_corner[:,1:-1]
                Yse_corner = Yse_corner[:,1:-1]
                Yne_corner = Yne_corner[:,1:-1]
                Ynw_corner = Ynw_corner[:,1:-1]
                Ysw_corner = Ysw_corner[:,1:-1]
        
                lon = XC.flatten()
                lat = YC.flatten()
                lon_corner1 = Xse_corner.flatten()
                lon_corner2 = Xsw_corner.flatten()
                lon_corner3 = Xnw_corner.flatten()
                lon_corner4 = Xne_corner.flatten()
                lat_corner1 = Yse_corner.flatten()
                lat_corner2 = Ysw_corner.flatten()
                lat_corner3 = Ynw_corner.flatten()
                lat_corner4 = Yne_corner.flatten()
            
            elif tile == 27 :
                # Remove edge cells without corner information
                XC = XC[0:-1,0:-1]
                YC = YC[0:-1,0:-1]
                Xse_corner = Xse_corner[0:-1,0:-1]
                Xne_corner = Xne_corner[0:-1,0:-1]
                Xnw_corner = Xnw_corner[0:-1,0:-1]
                Xsw_corner = Xsw_corner[0:-1,0:-1]
                Yse_corner = Yse_corner[0:-1,0:-1]
                Yne_corner = Yne_corner[0:-1,0:-1]
                Ynw_corner = Ynw_corner[0:-1,0:-1]
                Ysw_corner = Ysw_corner[0:-1,0:-1]

                ln = XC.flatten()
                lon = np.concatenate((lon, ln))
                lt = YC.flatten()
                lat = np.concatenate((lat,lt))
                lnc1 = Xse_corner.flatten()
                lon_corner1 = np.concatenate((lon_corner1,lnc1))
                lnc2 = Xsw_corner.flatten()
                lon_corner2 = np.concatenate((lon_corner2,lnc2))
                lnc3 = Xnw_corner.flatten()
                lon_corner3 = np.concatenate((lon_corner3,lnc3))
                lnc4 = Xne_corner.flatten()
                lon_corner4 = np.concatenate((lon_corner4,lnc4))
                
                ltc1 = Yse_corner.flatten()
                lat_corner1 = np.concatenate((lat_corner1,ltc1))
                ltc2 = Ysw_corner.flatten()
                lat_corner2 = np.concatenate((lat_corner2,ltc2))
                ltc3 = Ynw_corner.flatten()
                lat_corner3 = np.concatenate((lat_corner3,ltc3))
                ltc4 = Yne_corner.flatten()
                lat_corner4 = np.concatenate((lat_corner4,ltc4))

                # Create scrip file for combined unstructured grid
                grid_corner_lon = np.zeros((len(lon_corner1),4))
                grid_corner_lat = np.zeros((len(lat_corner1),4))
                grid_center_lon = np.zeros((len(lon),))
                grid_center_lat = np.zeros((len(lat),))
                grid_corner_lon[:,0] = lon_corner1
                grid_corner_lon[:,1] = lon_corner2
                grid_corner_lon[:,2] = lon_corner3
                grid_corner_lon[:,3] = lon_corner4
                grid_corner_lat[:,0] = lat_corner1
                grid_corner_lat[:,1] = lat_corner2
                grid_corner_lat[:,2] = lat_corner3
                grid_corner_lat[:,3] = lat_corner4
                grid_center_lon = lon
                grid_center_lat = lat
                grid_imask = np.ones((len(lon),)) # assume using all points
                grid_dims = np.array(lon.shape)

                scrip = xr.Dataset()
                scrip['grid_corner_lon'] = xr.DataArray(grid_corner_lon.astype('float64'),
                                                        dims=('grid_size','grid_corners'),
                                                        attrs={'units':'degrees'})
                scrip['grid_corner_lat'] = xr.DataArray(grid_corner_lat.astype('float64'),
                                                        dims=('grid_size','grid_corners'),
                                                        attrs={'units':'degrees'})
                scrip['grid_center_lon'] = xr.DataArray(grid_center_lon.astype('float64'),
                                                        dims=('grid_size'),
                                                        attrs={'units':'degrees'})
                scrip['grid_center_lat'] = xr.DataArray(grid_center_lat.astype('float64'),
                                                        dims=('grid_size'),
                                                        attrs={'units':'degrees'})
                scrip['grid_imask'] = xr.DataArray(grid_imask.astype('int32'),
                                                        dims=('grid_size'),
                                                        attrs={'units':'unitless'})
                scrip['grid_dims'] = xr.DataArray(grid_dims.astype('int32'),
                                                        dims=('grid_rank'),
                                                        attrs={'units':'unitless'})
                self.ecco_scripfile = 'ecco_combined_unstructured.scrip.nc'
                scrip.to_netcdf(self.ecco_scripfile)
                scrip.close()
                subprocess.run(["ncatted", "-a", "_FillValue,,d,,", self.ecco_scripfile])
                print("ECCO scrip file complete")

def main():
        run = eccoToMaliInterp()
        
        run.create_ECCO_scrip_file()
        
        run.merge_MITgcm_files_to_unstruct_grid()

        run.remap_ecco_to_mali()

if __name__ == "__main__":
    main()


