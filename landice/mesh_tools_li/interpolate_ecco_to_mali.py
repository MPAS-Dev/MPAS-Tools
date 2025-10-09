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
from pyproj import Proj, Transformer
from argparse import ArgumentParser
import time
import cftime
import glob

class eccoToMaliInterp:
    
    def __init__(self):
        print("Gathering Information ... ")
        parser = ArgumentParser(
                prog='interpolate_ecco_to_mali.py',
                description='Interpolates ECCO reanalysis data onto a MALI mesh')
        parser.add_argument("--eccoDir", dest="eccoDir", required=True, help="Directory where ECCO files are stored. Should be only netcdf files in that directory", metavar="FILENAME")
        parser.add_argument("--maliMesh", dest="maliMesh", required=True, help="MALI mesh to be interpolated onto", metavar="FILENAME")
        parser.add_argument("--tileNums", type=int, nargs='+', dest="tileNums", required=True, help="ECCO grid regional tile numbers to pull data from")
        parser.add_argument("--eccoVars", type=str, nargs='+', dest="eccoVars", required=True, help="ECCO variables to interpolate, current options are any combination of 'THETA', 'SALT', 'SIarea'")
        parser.add_argument("--maxDepth", type=int, dest="maxDepth", default=1500, help="Maximum depth in meters to include in the interpolation")
        parser.add_argument("--ntasks", dest="ntasks", type=str, help="Number of processors to use with ESMF_RegridWeightGen", default='128')
        parser.add_argument("--method", dest="method", type=str, help="Remapping method, either 'bilinear' or 'neareststod'", default='conserve')
        parser.add_argument("--outFile", dest="outputFile", help="Desired name of output file", metavar="FILENAME", default="ecco_to_mali_remapped.nc")
        parser.add_argument("--meshVars", dest="includeMeshVars", help="whether to include mesh variables in resulting MALI forcing file", action='store_true')
        args = parser.parse_args()
        self.options = args

        self.mapping_file_name = f'ecco_to_mali_mapping_{self.options.method}.nc'

    def merge_MITgcm_files_to_unstruct_grid(self):
        ncfiles = glob.glob(os.path.join(self.options.eccoDir, "*.nc"))
        saltBool = 0
        thetaBool = 0
        siaBool = 0
    
        #loop through MITgcm tile files and create unstructured arrays for each variable in eccoVars
        for tile in self.options.tileNums:
            print(f"Processing Tile {tile}")
            o3dmBool = 0
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
                            
                        sal = ds_ecco['SALT'][:].values
                        
                        # ECCO 'land' variable is equivalent to MALI's orig3dOceanMask, despite the counterintuitive name.
                        # Only need to call this once per tile, but needs to be using a theta/salinity file because it includes depth 
                        # No need to call this if only interpolating SIarea
                        if o3dmBool == 0:
                            orig3dOceanMask = ds_ecco['land'].values
                            od3mBool = 1
                        
                        # set MALI invalid value
                        boolZ,boolX,boolY = np.where(orig3dOceanMask == 0)
                        for iZ,iX,iY in zip(boolZ,boolX,boolY):
                            sal[:,iZ,iX,iY] = 1e36
                        
                        saltBool = 1       

                    # Repeat for potential temperature if needed
                    elif 'THETA' in file and 'THETA' in self.options.eccoVars:       
                        validFile = self.options.eccoDir + '/THETA.' + paddedTile + '.nc'
                        ds_ecco = xr.open_dataset(validFile)
                        
                        if saltBool == 0 and thetaBool == 0:
                            z = ds_ecco['dep'].values
                            
                        theta = ds_ecco['THETA'][:].values

                        # ECCO 'land' variable is equivalent to MALI's orig3dOceanMask, despite the counterintuitive name.
                        # Only need to call this once per tile, but needs to be using a theta/salinity file because it includes depth 
                        # No need to call this if only interpolating SIarea
                        if o3dmBool == 0:
                            orig3dOceanMask = ds_ecco['land'].values
                            od3mBool = 1
                        
                        # set MALI invalid value
                        boolZ,boolX,boolY = np.where(orig3dOceanMask == 0)
                        for iZ,iX,iY in zip(boolZ,boolX,boolY):
                            theta[:,iZ,iX,iY] = 1e36
                        
                        thetaBool = 1
                        
                    # Repeat for sea ice fraction if needed
                    elif 'SIarea' in file and 'SIarea' in self.options.eccoVars:
                        validFile = self.options.eccoDir + '/SIarea.' + paddedTile + '.nc'
                        ds_ecco = xr.open_dataset(validFile)
                        
                        sia = ds_ecco['SIarea'][:].values
                    
                        siaBool = 1
                    else :
                        raise ValueError(f"eccoDir: {self.options.eccoDir}, paddedTile: {paddedTile}")

            if ('SALT' in self.options.eccoVars and saltBool == 0):
                raise ValueError(f"SALT netcdf file for tile {tile} not found in {self.options.eccoDir}")

            if ('THETA' in self.options.eccoVars and thetaBool == 0):
                raise ValueError(f"THETA netcdf file for tile {tile} not found in {self.options.eccoDir}")
                
            if ('SIarea' in self.options.eccoVars and siaBool == 0):
                raise ValueError(f"SIarea netcdf file for tile {tile} not found in {self.options.eccoDir}")
                
                            
            # Doesn't matter which file these come from, just that there is only one set per tile, so just use most
            # recent valid file in each tile
        
            lon = ds_ecco['lon'][:].values.flatten()
            lat = ds_ecco['lat'][:].values.flatten()

            # combine into new netcdf. Need to load depth and time from original ecco files. Assuming time/depth is consistent
            # across all files, just load info from most recent
            
            # Establish dataset, and fill in standard dimensions. other dimensions are added on an as-needed basis
            ds_tile = xr.Dataset()
            ds_tile.expand_dims(["nCells", "Time"])
            
            # Translate ECCO time to MALI xtime format as save to dataset
            time = ds_ecco['tim'].values
            xtime = [pd.to_datetime(dt).strftime("%Y-%m-%d_%H:%M:%S").ljust(64) for dt in time]
            
            DA_xtime = xr.DataArray(np.array(xtime).astype('S64'), dims=("Time"))
            DA_xtime.encoding.update({"char_dim_name":"StrLen"})
            DA_xtime.attrs['long_name'] = "model time, with format 'YYYY-MM-DD_HH:MM:SS'"
            DA_xtime.attrs['standard_name'] = 'time'
            ds_tile['xtime'] = DA_xtime
            
            # save unstructured variables with MALI/MPAS-O names 
            if 'z' in locals():
                ds_tile.expand_dims("nISMIP6OceanLayers")
                ds_tile.expand_dims("TWO")
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
                ds_tile['ismip6shelfMelt_zOcean'] = DA_z

                DA_zBnds = xr.DataArray(-zBnds.astype('float64'), dims=("nISMIP6OceanLayers", "TWO")) #MALI wants negative depth coordinates
                DA_zBnds.attrs['long_name'] = 'Bounds for ismip6 ocean layers'
                DA_zBnds.attrs['units'] = 'm'
                ds_tile['ismip6shelfMelt_zBndsOcean'] = DA_zBnds
            
            if 'theta' in locals():
                print(f"z shape: {z.shape}")
                print(f"theta shape: {theta.shape}")
                theta = theta[:,0:len(z),:,:]
                print(f"new theta shape: {theta.shape}")
                DA_theta = xr.DataArray(theta.astype('float64'), dims=("Time","nISMIP6OceanLayers", "X", "Y"))
                DA_theta.attrs['long_name'] = 'Potential_Temperature'
                DA_theta.attrs['units'] = 'degC'
                ds_tile['oceanTemperature'] = DA_theta
                    
            if 'sal' in locals():
                sal = sal[:,0:len(z),:,:]
                DA_salt = xr.DataArray(sal.astype('float64'), dims=("Time","nISMIP6OceanLayers", "X", "Y"))
                DA_salt.attrs['long_name'] = 'Salinity'
                DA_salt.attrs['units'] = 'psu'
                ds_tile['oceanSalinity'] = DA_salt
                    
            if 'sia' in locals():
                DA_sia = xr.DataArray(sia.astype('float64'), dims=("Time", "X", "Y"))
                DA_sia.attrs['long_name'] = 'Sea ice fractional ice-covered area [0 to 1]'
                DA_sia.attrs['units'] = 'm2/m2'
                ds_tile['iceAreaCell'] = DA_sia

            if 'orig3dOceanMask' in locals():
                orig3dOceanMask = orig3dOceanMask[0:len(z),:]
                DA_o3dm = xr.DataArray(orig3dOceanMask.astype('int32'), dims=("nISMIP6OceanLayers", "X", "Y"))
                DA_o3dm.attrs['long_name'] = ("3D mask of original valid ocean data.  Because it is 3d," 
                                             " it can include the ocean domain inside ice-shelf cavities."
                                             " orig3dOceanMask is altered online to include ice-dammed inland seas")
                ds_tile['orig3dOceanMask'] = DA_o3dm 

            DA_lon = xr.DataArray(lon.astype('float64'), dims=("X"))
            DA_lon.attrs['long_name'] = 'longitude'
            DA_lon.attrs['standard_name'] = 'longitude'
            DA_lon.attrs['units'] = 'degree_east'
            ds_tile['lon'] = DA_lon # Keep original lon/lat names of coordinates for projection transformer

            DA_lat = xr.DataArray(lat.astype('float64'), dims=("Y"))
            DA_lat.attrs['long_name'] = 'latitude'
            DA_lat.attrs['standard_name'] = 'latitude'
            DA_lat.attrs['units'] = 'degree_north'
            ds_tile['lat'] = DA_lat
            
            #convert ecco mesh to polar stereographic
            #use gis-gimp projection, defined in mpas_tools.landice.projections
            proj_string = (
            "+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 "
            "+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84"
            )
            transformer = Transformer.from_crs("EPSG:4326", proj_string, always_xy=True)
            x,y = transformer.transform(lon,lat)
            
            DA_x['x'] = xr.DataArray(x.astype('float64'), dims=("X"))
            ds_tile['x'] = DA_x
            
            DA_y['y'] = xr.DataArray(y.astype('float64'), dims=("Y"))
            ds_tile['y'] = DA_y

            ecco_tile = f"ecco_combined_tile{tile}.nc"
            ds_tile.to_netcdf(ecco_tile)
            ds_tile.close()

            print("Creating Scrip Files")
            filename = os.path.basename(ecco_tile)
            base, ext = os.path.splitext(filename)
            ecco_scripfile = base + '.scrip' + ext
            
            filename = os.path.basename(self.options.maliMesh)
            base, ext = os.path.splitext(filename)
            mali_scripfile = base + '.scrip' + ext
            
            args = ['create_scrip_file_from_planar_rectangular_grid',
                '-i', ecco_tile,
                '-s', ecco_scripfile,
                '-p', 'gis-gimp',
                '-r', '2']
            check_call(args)
            
            scrip_from_mpas(self.options.maliMesh, mali_scripfile)

            #create mapping file
            print("Creating Mapping File")
            args = ['srun',
                    '-n', self.options.ntasks, 'ESMF_RegridWeightGen',
                    '--source', ecco_scripfile,
                    '--destination', mali_scripfile,
                    '--weight', self.mapping_file_name,
                    '--method', self.options.method,
                    '--netcdf4',
                    "--dst_regional", "--src_regional", '--ignore_unmapped']
            check_call(args)
    
            #remap
            print("Start remapping ... ")

            subprocess.run(["ncatted", "-a", "_FillValue,,d,,", ecco_tile])

            print("ncremap ...")
            st = time.time()
            # remap the input data
            args = ["ncremap",
                    "-i", ecco_tile,
                    "-o", self.options.outputFile,
                    "-m", self.mapping_file_name]
            check_call(args)

            subprocess.run(["ncpdq", "-O", "-a", "Time,nCells,nISMIP6OceanLayers", self.options.outputFile, self.options.outputFile])
            nd = time.time()
            tm = nd - st
            print(f"ncremap completed: {tm} seconds")
            print("Remapping completed")

        # Transfer MALI mesh variables if needed
        if self.options.includeMeshVars:
            ds_mali = xr.open_dataset(self.options.maliMesh, decode_times=False, decode_cf=False)
            ds_out = xr.open_dataset(self.options.outputFile, decode_times=False, decode_cf=False)
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
            ds_out.to_netcdf(self.options.outputFile, mode='w')
            ds_out.close()
            ds_mali.close()
        
def main():
        run = eccoToMaliInterp()
        
        run.merge_MITgcm_files_to_unstruct_grid()

if __name__ == "__main__":
    main()


