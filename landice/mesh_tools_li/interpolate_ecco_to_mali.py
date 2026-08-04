#!/usr/bin/env python
'''
This script interpolates ECCO ASTE R1 reanalysis data
(https://arcticdata.io/data/10.18739/A2CV4BS5K/) onto a MALI mesh.

Ocean temperature and salinity are automatically interpolated and the user has
the option to also interpolate sea ice fraction. Interpolation script will
generate MALI input field icebergFjordMask from sea ice area if using this
option.

ECCO grid files are unconventional for MITgcm simulations and do not contain
adequate grid cell corner information for all tile edge cells (necessary for
creating scrip file and remapping). This script includes seam-aware stitching
for the Greenland ASTE subset (tiles 5, 11, 12, 14, 15, 27).
'''

import os
import shutil
import xarray as xr
import subprocess
import numpy as np
import pandas as pd
from mpas_tools.logging import check_call
from mpas_tools.scrip.from_mpas import scrip_from_mpas
from argparse import ArgumentParser
import time
import cftime
import json
try:
    from shapely.geometry import shape
    from shapely.prepared import prep
    from shapely.geometry import Point
except ImportError:
    shape = None
    prep = None
    Point = None

class eccoToMaliInterp:
    
    def __init__(self):
        print("Gathering Information ... ")
        parser = ArgumentParser(
                prog='interpolate_ecco_to_mali.py',
                description=__doc__)
        parser.add_argument("--eccoDir", dest="eccoDir", required=True, help="Directory where ECCO files are stored. Should be only netcdf files in that directory", metavar="FILENAME")
        parser.add_argument("--maliMesh", dest="maliMesh", required=True, help="MALI mesh to be interpolated onto", metavar="FILENAME")
        parser.add_argument("--sia", dest="sia", help="Option to remap sea ice area in addition to temp/salt. Needed for icebergFjordMask", action='store_true')
        parser.add_argument("--geojson", dest="geojson", help="Optional geojson used to permanently define icebergFjordMask within a region")
        parser.add_argument("--maxDepth", type=int, dest="maxDepth", default=1500, help="Maximum depth in meters to include in the interpolation")
        parser.add_argument("--ntasks", dest="ntasks", type=str, help="Number of processors to use with ESMF_RegridWeightGen", default='128')
        parser.add_argument("--method", dest="method", type=str, help="Remapping method, either 'bilinear' or 'neareststod'", default='bilinear')
        parser.add_argument("--outFile", dest="outputFile", help="Desired name of output file", metavar="FILENAME", default="ecco_to_mali_remapped.nc")
        parser.add_argument("--meshVars", dest="includeMeshVars", help="whether to include mesh variables in resulting MALI forcing file", action='store_true')
        parser.add_argument("--iceAreaThresh", dest="iceAreaThresh", type=float, help="Threshold sea ice fraction used to define icebergFjordMask.", default=0.5)
        parser.add_argument("--extrapIters", dest="extrapIters", type=str, help="Number of sea ice extrapolation iterations when regridding", default='175')
        parser.add_argument("--keepTmpFiles", dest="keepTmpFiles", help="Whether to keep temporary netcdf files created throughout script. Useful for debugging", action="store_true")
        args = parser.parse_args()
        self.options = args

        self.mapping_file_name = f'ecco_to_mali_mapping_{self.options.method}.nc'
        self.supportedTiles = {5, 11, 12, 14, 15, 27}
        self.tileNums = sorted(self.supportedTiles)
        # Hardcoded seams for Greenland ASTE subset: tile_edge -> neighbor_edge.
        self.greenlandSeams = {
            5: {'right': (14, 'bottom', True)},
            11: {'right': (14, 'left', False), 'top': (12, 'bottom', False)},
            12: {'right': (15, 'left', False)},
            14: {'right': (27, 'bottom', True), 'top': (15, 'bottom', False)},
            27: {'right': (5, 'bottom', True)}
        }
        self._tileGeom = None

    @staticmethod
    def _edge_vector(xg, yg, edge_name):
        if edge_name == 'left':
            return xg[0, :].copy(), yg[0, :].copy()
        if edge_name == 'right':
            return xg[-1, :].copy(), yg[-1, :].copy()
        if edge_name == 'bottom':
            return xg[:, 0].copy(), yg[:, 0].copy()
        if edge_name == 'top':
            return xg[:, -1].copy(), yg[:, -1].copy()
        raise ValueError(f"Unknown edge name: {edge_name}")

    def _build_greenland_tile_geometry(self):
        if self._tileGeom is not None:
            return self._tileGeom

        grids = {}
        for tile in self.tileNums:
            padded_tile = str(tile).zfill(4)
            grid_file = os.path.join(self.options.eccoDir, f'GRID.{padded_tile}.nc')
            if not os.path.exists(grid_file):
                raise ValueError(f"GRID netcdf file for tile {tile} not found in {self.options.eccoDir}")
            with xr.open_dataset(grid_file) as grd:
                grids[tile] = {
                    'XC': grd['XC'][:].values,
                    'YC': grd['YC'][:].values,
                    'XG': grd['XG'][:].values,
                    'YG': grd['YG'][:].values
                }

        tile_geom = {}
        for tile in self.tileNums:
            xc = grids[tile]['XC']
            yc = grids[tile]['YC']
            xg = grids[tile]['XG']
            yg = grids[tile]['YG']
            nx, ny = xc.shape

            right_vec_lon = None
            right_vec_lat = None
            top_vec_lon = None
            top_vec_lat = None

            tile_seams = self.greenlandSeams.get(tile, {})

            if 'right' in tile_seams:
                neighbor_tile, neighbor_edge, reverse_order = tile_seams['right']
                if neighbor_tile in grids:
                    rv_lon, rv_lat = self._edge_vector(grids[neighbor_tile]['XG'], grids[neighbor_tile]['YG'], neighbor_edge)
                    if reverse_order:
                        rv_lon = rv_lon[::-1]
                        rv_lat = rv_lat[::-1]
                    right_vec_lon, right_vec_lat = rv_lon, rv_lat

            if 'top' in tile_seams:
                neighbor_tile, neighbor_edge, reverse_order = tile_seams['top']
                if neighbor_tile in grids:
                    tv_lon, tv_lat = self._edge_vector(grids[neighbor_tile]['XG'], grids[neighbor_tile]['YG'], neighbor_edge)
                    if reverse_order:
                        tv_lon = tv_lon[::-1]
                        tv_lat = tv_lat[::-1]
                    top_vec_lon, top_vec_lat = tv_lon, tv_lat

            i_cells = []
            j_cells = []
            lon_cells = []
            lat_cells = []
            c1_lon = []
            c2_lon = []
            c3_lon = []
            c4_lon = []
            c1_lat = []
            c2_lat = []
            c3_lat = []
            c4_lat = []

            for i in range(nx):
                for j in range(ny):
                    se_lon = xg[i, j]
                    se_lat = yg[i, j]

                    if i + 1 < nx:
                        sw_lon = xg[i + 1, j]
                        sw_lat = yg[i + 1, j]
                    elif right_vec_lon is not None:
                        sw_lon = right_vec_lon[j]
                        sw_lat = right_vec_lat[j]
                    else:
                        continue

                    if j + 1 < ny:
                        ne_lon = xg[i, j + 1]
                        ne_lat = yg[i, j + 1]
                    elif top_vec_lon is not None:
                        ne_lon = top_vec_lon[i]
                        ne_lat = top_vec_lat[i]
                    else:
                        continue

                    if i + 1 < nx and j + 1 < ny:
                        nw_lon = xg[i + 1, j + 1]
                        nw_lat = yg[i + 1, j + 1]
                    elif i + 1 == nx and j + 1 < ny and right_vec_lon is not None:
                        nw_lon = right_vec_lon[j + 1]
                        nw_lat = right_vec_lat[j + 1]
                    elif i + 1 < nx and j + 1 == ny and top_vec_lon is not None:
                        nw_lon = top_vec_lon[i + 1]
                        nw_lat = top_vec_lat[i + 1]
                    else:
                        # Corner cell where both +i and +j need cross-face corner.
                        continue

                    i_cells.append(i)
                    j_cells.append(j)
                    lon_cells.append(xc[i, j])
                    lat_cells.append(yc[i, j])
                    c1_lon.append(se_lon)
                    c2_lon.append(sw_lon)
                    c3_lon.append(nw_lon)
                    c4_lon.append(ne_lon)
                    c1_lat.append(se_lat)
                    c2_lat.append(sw_lat)
                    c3_lat.append(nw_lat)
                    c4_lat.append(ne_lat)

            tile_geom[tile] = {
                'i': np.array(i_cells, dtype=np.int32),
                'j': np.array(j_cells, dtype=np.int32),
                'center_lon': np.array(lon_cells, dtype=np.float64),
                'center_lat': np.array(lat_cells, dtype=np.float64),
                'corner_lon': np.stack((c1_lon, c2_lon, c3_lon, c4_lon), axis=1).astype(np.float64),
                'corner_lat': np.stack((c1_lat, c2_lat, c3_lat, c4_lat), axis=1).astype(np.float64)
            }

        self._tileGeom = tile_geom
        return tile_geom

    def merge_MITgcm_files_to_unstruct_grid(self):
        tile_geom = self._build_greenland_tile_geometry()

        sal_tiles = []
        theta_tiles = []
        sia_tiles = []
        mask_tiles = []
        z = None
        time_vals = None

        # Loop through selected MITgcm tile files and create unstructured arrays.
        for tile in self.tileNums:
            paddedTile = str(tile).zfill(4)
            print(f"Processing tile {paddedTile}")

            salt_file = os.path.join(self.options.eccoDir, f"SALT.{paddedTile}.nc")
            theta_file = os.path.join(self.options.eccoDir, f"THETA.{paddedTile}.nc")
            sia_file = os.path.join(self.options.eccoDir, f"SIarea.{paddedTile}.nc")

            if not os.path.exists(salt_file):
                raise ValueError(f"SALT netcdf file for tile {tile} not found in {self.options.eccoDir}")
            if not os.path.exists(theta_file):
                raise ValueError(f"THETA netcdf file for tile {tile} not found in {self.options.eccoDir}")
            if self.options.sia and not os.path.exists(sia_file):
                raise ValueError(f"SIarea netcdf file for tile {tile} not found in {self.options.eccoDir}")

            with xr.open_dataset(salt_file) as ds_salt:
                if z is None:
                    z = ds_salt['dep'].values
                if time_vals is None:
                    time_vals = ds_salt['tim'].values
                i_idx = tile_geom[tile]['i']
                j_idx = tile_geom[tile]['j']
                sal_tile = ds_salt['SALT'][:].values
                nt, nz, nx, ny = np.shape(sal_tile)
                sal_tiles.append(sal_tile[:, :, i_idx, j_idx])

                # ECCO 'land' is equivalent to MALI's orig3dOceanMask.
                mask_tile = ds_salt['land'][:].values
                mask_tiles.append(mask_tile[:, i_idx, j_idx])

            with xr.open_dataset(theta_file) as ds_theta:
                theta_tile = ds_theta['THETA'][:].values
                theta_tiles.append(theta_tile[:, :, i_idx, j_idx])

            if self.options.sia:
                with xr.open_dataset(sia_file) as ds_sia_in:
                    sia_tile = ds_sia_in['SIarea'][:].values
                    sia_tiles.append(sia_tile[:, i_idx, j_idx])

        sal = np.concatenate(sal_tiles, axis=2)
        theta = np.concatenate(theta_tiles, axis=2)
        orig3dOceanMask = np.concatenate(mask_tiles, axis=1)
        nt = sal.shape[0]
        if self.options.sia:
            sia = np.concatenate(sia_tiles, axis=1)

        # Set MALI invalid value over non-ocean cells.
        ocean_mask = orig3dOceanMask[np.newaxis, :, :] > 0
        sal = np.where(ocean_mask, sal, 1e36)
        theta = np.where(ocean_mask, theta, 1e36)

        # combine into new netcdf. Need to load depth and time from original ecco files. Assuming time/depth is consistent
        # across all files, just load info from most recent
        
        # Establish dataset, and fill in standard dimensions. other dimensions are added on an as-needed basis
        ds_unstruct = xr.Dataset()
        ds_unstruct.expand_dims(["nCells", "Time"])
        
        # Save ECCO time variable as datetime string. Will be converted into xtime format when saving final output file
        # Subtract one month from each time step. Ecco to output every month as an average of the previous month. We
        # want to force MALI with ECCO data from the same month (avoid 1 month lag)
        time = time_vals
        self.xtime_str = [(pd.to_datetime(dt) - pd.DateOffset(months=1)).strftime("%Y-%m-%d_%H:%M:%S").ljust(64) for dt in time]
    
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

        if self.options.sia:
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
            ds_sia = xr.open_dataset(siaRemappedOutFile)
            iac = ds_sia['iceAreaCell'].values
            icebergFjordMask = np.zeros(iac.shape)
            # Find cells where sea ice fraction exceeds threshold
            seaice_mask = iac > self.options.iceAreaThresh 
            icebergFjordMask[seaice_mask] = 1

            # Use geojson to add any additional cells to icebergFjordMask
            if self.options.geojson is not None:
                if shape is None or prep is None or Point is None:
                    raise ImportError("shapely is required when using --geojson")
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
        ds_out = xr.open_dataset(remappedOutFile, decode_times=False, decode_cf=False, mode='r+')
        o3dm = ds_out['orig3dOceanMask']
        temp = ds_out['oceanTemperature']
        sal = ds_out['oceanSalinity']

        # Account for rounding error. Good cells are not exactly 1 and temp/sal invalid values may be averaged out.
        # It is otherwise possible to get some invalid temp/sal values align with good orig3dOceanMask values after remapping
        mask = (o3dm > 0.9999) & (temp < 100) & (sal < 100)
        ds_out['orig3dOceanMask'] = xr.DataArray(mask.astype('int32'), dims=("Time", "nCells", "nISMIP6OceanLayers"))
        ds_out['oceanTemperature'] = temp.where(mask, 1e36).astype('float64')
        ds_out['oceanSalinity'] = sal.where(mask, 1e36).astype('float64')
        if self.options.sia:
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
            ds_out.to_netcdf(self.options.outputFile, unlimited_dims=['Time'])
            subprocess.run(["ncatted", "-a", "_FillValue,,d,,", self.options.outputFile])

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
            ds_out.to_netcdf(self.options.outputFile, unlimited_dims=['Time'])
            ds_out.close()
            ds_mali.close()
            subprocess.run(["ncatted", "-a", "_FillValue,,d,,", self.options.outputFile])
        
        # Remove temporary output files
        if not self.options.keepTmpFiles :
            for filename in os.listdir("."):
                if "tmp" in filename and os.path.isfile(filename):
                    os.remove(filename)

        # Convert to netcdf3 64-bit offset. Much faster to do this with ncks than create the netcdf3 originally
        # with xarray.to_netcdf. 
        nc64offset = self.options.outputFile[0:-3] + '.64bitOffset.nc'
        subprocess.run(["ncks", "-O", "-6", self.options.outputFile, nc64offset])

    def create_ECCO_scrip_file(self):
        
        print("Creating ECCO scrip file")
        tile_geom = self._build_greenland_tile_geometry()

        lon_parts = []
        lat_parts = []
        lon_corner1_parts = []
        lon_corner2_parts = []
        lon_corner3_parts = []
        lon_corner4_parts = []
        lat_corner1_parts = []
        lat_corner2_parts = []
        lat_corner3_parts = []
        lat_corner4_parts = []

        for tile in self.tileNums:
            lon_parts.append(tile_geom[tile]['center_lon'])
            lat_parts.append(tile_geom[tile]['center_lat'])
            lon_corner1_parts.append(tile_geom[tile]['corner_lon'][:, 0])
            lon_corner2_parts.append(tile_geom[tile]['corner_lon'][:, 1])
            lon_corner3_parts.append(tile_geom[tile]['corner_lon'][:, 2])
            lon_corner4_parts.append(tile_geom[tile]['corner_lon'][:, 3])
            lat_corner1_parts.append(tile_geom[tile]['corner_lat'][:, 0])
            lat_corner2_parts.append(tile_geom[tile]['corner_lat'][:, 1])
            lat_corner3_parts.append(tile_geom[tile]['corner_lat'][:, 2])
            lat_corner4_parts.append(tile_geom[tile]['corner_lat'][:, 3])

        lon = np.concatenate(lon_parts)
        lat = np.concatenate(lat_parts)
        lon_corner1 = np.concatenate(lon_corner1_parts)
        lon_corner2 = np.concatenate(lon_corner2_parts)
        lon_corner3 = np.concatenate(lon_corner3_parts)
        lon_corner4 = np.concatenate(lon_corner4_parts)
        lat_corner1 = np.concatenate(lat_corner1_parts)
        lat_corner2 = np.concatenate(lat_corner2_parts)
        lat_corner3 = np.concatenate(lat_corner3_parts)
        lat_corner4 = np.concatenate(lat_corner4_parts)

        # Create scrip file for combined unstructured grid.
        grid_corner_lon = np.zeros((len(lon_corner1), 4))
        grid_corner_lat = np.zeros((len(lat_corner1), 4))
        grid_corner_lon[:, 0] = lon_corner1
        grid_corner_lon[:, 1] = lon_corner2
        grid_corner_lon[:, 2] = lon_corner3
        grid_corner_lon[:, 3] = lon_corner4
        grid_corner_lat[:, 0] = lat_corner1
        grid_corner_lat[:, 1] = lat_corner2
        grid_corner_lat[:, 2] = lat_corner3
        grid_corner_lat[:, 3] = lat_corner4
        grid_imask = np.ones((len(lon),))
        grid_dims = np.array(lon.shape)

        scrip = xr.Dataset()
        scrip['grid_corner_lon'] = xr.DataArray(grid_corner_lon.astype('float64'),
                                                dims=('grid_size', 'grid_corners'),
                                                attrs={'units': 'degrees'})
        scrip['grid_corner_lat'] = xr.DataArray(grid_corner_lat.astype('float64'),
                                                dims=('grid_size', 'grid_corners'),
                                                attrs={'units': 'degrees'})
        scrip['grid_center_lon'] = xr.DataArray(lon.astype('float64'),
                                                dims=('grid_size'),
                                                attrs={'units': 'degrees'})
        scrip['grid_center_lat'] = xr.DataArray(lat.astype('float64'),
                                                dims=('grid_size'),
                                                attrs={'units': 'degrees'})
        scrip['grid_imask'] = xr.DataArray(grid_imask.astype('int32'),
                                           dims=('grid_size'),
                                           attrs={'units': 'unitless'})
        scrip['grid_dims'] = xr.DataArray(grid_dims.astype('int32'),
                                          dims=('grid_rank'),
                                          attrs={'units': 'unitless'})
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


