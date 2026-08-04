#!/usr/bin/env python

from pathlib import Path
import xarray as xr
import numpy as np
from argparse import ArgumentParser
from datetime import datetime
from pyproj import Proj, Transformer
import json
from shapely.geometry import shape, Point
from scipy.spatial import cKDTree
import cftime
from time import perf_counter
import pandas as pd
from matplotlib import pyplot as plt

class validateMaliExtrap:
    def __init__(self):
        parser = ArgumentParser(
                prog = 'validate_mali_extrap.py',
                description = 'Synthesizes a CTD datasets in the Uummannaq/Disko region and calculates a Willmot Skill Score for a given MALI thermal forcing paramterization.')
        parser.add_argument("--maliDiag", dest='maliDiag', required=True, help="MALI output file to validate")
        parser.add_argument("--maliForcing", dest='maliForcing', required=True, help="MALI ocean forcing file")
        parser.add_argument("--HadleyDir", dest="HadleyDir", help="Directory with Hadley EN4 Casts")
        parser.add_argument("--geojson", dest="geojson", help="geojson defining bounds of domain (lat/lon)")
        args = parser.parse_args()
        self.options = args

    def process_Hadley_EN4(self):

        files = list(Path(self.options.HadleyDir).rglob("*.nc"))

        depth = []
        thermal_forcing = []
        time = []
        lat = []
        lon = []

        t0 = perf_counter()
        ct = 0
        processedFiles = np.arange(1,len(files),100)

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(1,1,1)
        
        for fpath in files:
            ct = ct + 1
            if np.any(processedFiles == ct):
                print(f"processed {ct} of {len(files)} files")

            ds = xr.open_dataset(fpath, decode_cf=False, decode_times=False)

            z = ds["depth"].values
            depth.append(z)

            temp = ds["potential_temperature"].values
            salt = ds["practical_salinity"].values

            tf = _calc_thermal_forcing(temp, salt, -z)
            thermal_forcing.append(tf)

            attrs = ds.attrs
            time.append(datetime(
                int(attrs["year"]),
                int(attrs["month"]),
                int(attrs["day"]),
                int(attrs.get("hour", 0)),
                int(attrs.get("minute", 0))
            ))
            lat.append(attrs["latitude"])
            lon.append(attrs["longitude"])

            ds.close()

            ax.plot(tf,-z)
        
        fig.savefig('initialData.png')
        print(f"opening files took: {perf_counter() - t0:.1f} s")

        # Convert to (depth,profile) format and interpolate to MALI z coordinates
        ds_maliForcing = xr.open_dataset(self.options.maliForcing, decode_times=False, decode_cf=False)
        z_mali = np.abs(ds_maliForcing['ismip6shelfMelt_zOcean'][:].values)
        
        TF = np.full((len(thermal_forcing), z_mali.size), np.nan)

        print("interpolating ...")

        t0 = perf_counter()
        for i, (z, tf) in enumerate(zip(depth, thermal_forcing)):
            
            mask = np.isfinite(z) & np.isfinite(tf)
            if mask.sum() < 2:
                continue

            TF[i, :] = np.interp(
                z_mali,
                z[mask],
                tf[mask],
                left=np.nan,
                right=np.nan,
            )
        if np.all(np.isnan(TF)):
            print(f"TF is all NaNs")

        print(f"Interpolation took: {perf_counter() - t0:.1f} s")
        print(f"TF shape: {TF.shape}")
        en4lon = np.array(lon)
        en4lat = np.array(lat)
        time = pd.to_datetime(time) 
        print(f"lon shape: {en4lon.shape}")
        print(f"lat shape: {en4lat.shape}") 

        # Filter casts to exclude anything outside of MALI mesh. Needs geojson
        with open(self.options.geojson) as f:
            geojson = json.load(f)
        region = shape(geojson["features"][0]["geometry"])

        mask = np.array([
            region.contains(Point(lon, lat))
            for lon, lat in zip(en4lon, en4lat)
        ])

        en4lon = en4lon[mask]
        en4lat = en4lat[mask]
        TF = TF[mask,:]
        time = time[mask]

        # filter by time. Only want years overlapping with ECCO
        mask = (time.year >= 2002) & (time.year <= 2018)
        en4lon = en4lon[mask]
        en4lat = en4lat[mask]
        TF = TF[mask,:]
        time = time[mask]
        
        # coordinate conversion onto polar stereographic
        transformer = Transformer.from_crs("EPSG:4326", 
                                           "+proj=stere +lat_ts=70 +lat_0=90 +lon_0=315 +k_0=1 +x_0=0 +y_0=0 +ellps=WGS84", 
                                           always_xy=True)
        en4x, en4y = transformer.transform(en4lon, en4lat)
    
        ds_maliDiag = xr.open_dataset(self.options.maliDiag, decode_times=False, decode_cf=False)
        mali_xy = np.column_stack((ds_maliDiag['xCell'][:].values, ds_maliDiag['yCell'][:].values))
        tree = cKDTree(mali_xy)
        dist, ind = tree.query(np.column_stack((en4x, en4y)), k=1)

        # Mask ocean cells with original ECCO data vs extrapolated values
        o3dm = ds_maliForcing['orig3dOceanMask'][0, :, 0].values
        z = ds_maliDiag['bedTopography'][0, :].values
        oceanDataMask = np.zeros(ind.shape, dtype='int32')
        oceanDataMask[(z[ind] < 0) & (o3dm[ind] == 0)] = 1 # extrapolated
        oceanDataMask[(z[ind] < 0) & (o3dm[ind] == 1)] = 2 # original ecco data


        # prepare for WSS comparison
        tf_mali = ds_maliDiag['ismip6shelfMelt_3dThermalForcing'][:].values
        tf_mali[tf_mali > 100] = np.nan

        print("xtime")
        xtime = xr.DataArray(
            [cftime.datetime.strptime(x.tobytes().decode('utf-8').strip("\x00"), "%Y-%m-%d_%H:%M:%S")
             for x in ds_maliDiag.xtime.values],
            dims=("Time",)
        )

        # save variables
        self.yearMonth_mali = np.array((xtime.dt.year * 100 + xtime.dt.month).data).astype('int32')
        self.yearMonth_casts = np.array(time.year * 100 + time.month).astype('int32')
        self.corr_mali_ind = ind
        self.oceanDataMask = oceanDataMask
        self.tf_casts = TF
        self.tf_mali = tf_mali
        self.z = z_mali

        ds_maliDiag.close()
        ds_maliForcing.close()

    def calc_willmott_skill_score(self):
        print("calc_willmott")
        print(f"yearMonth_mali: {self.yearMonth_mali.shape}")
        print(f"yearMonth_ecco: {self.yearMonth_casts.shape}")
        print(f"corr_mali_ind: {self.corr_mali_ind.shape}")
        print(f"oceanDataMask: {self.oceanDataMask.shape}")
        print(f"tf_casts: {self.tf_casts.shape}")
        print(f"tf_mali: {self.tf_mali.shape}")
        print(f"z: {self.z.shape}")

        origCastMask = self.oceanDataMask == 2
        extrapCastMask = self.oceanDataMask == 1

        tf_o_c = self.tf_casts[origCastMask, :]
        tf_e_c = self.tf_casts[extrapCastMask, :]
        print(f"tf_o_c: {tf_o_c.shape}")
        print(f"tf_e_c: {tf_e_c.shape}")

        yearMonth_o_c = self.yearMonth_casts[origCastMask]
        yearMonth_e_c = self.yearMonth_casts[extrapCastMask]
        print(f"yearMonth_o_c: {yearMonth_o_c.shape}")
        print(f"yearMonth_e_c: {yearMonth_e_c.shape}")

        indMali_o_c = self.corr_mali_ind[origCastMask]
        indMali_e_c = self.corr_mali_ind[extrapCastMask]
        print(f"indMali_o_c: {indMali_o_c.shape}")
        print(f"indMali_e_c: {indMali_e_c.shape}")

        # average any casts with same yearMonth and indMali pairings
        def _average_ctd_profiles(yearMonth,indMali,tf):    
            pairs = np.column_stack((yearMonth, indMali))
            unique_pairs, inverse = np.unique(pairs, axis=0, return_inverse=True)
            
            K = len(unique_pairs)
            Z = tf.shape[1]
            tf_avg = np.zeros((K, Z))
            
            counts = np.bincount(inverse)
            
            np.add.at(tf_avg, inverse, tf)
            tf_avg /= counts[:, None]

            yearMonth_condensed = unique_pairs[:, 0]
            indMali_condensed = unique_pairs[:, 1]
            
            return yearMonth_condensed, indMali_condensed, tf_avg

        yearMonth_o_c, indMali_o_c, tf_o_c = _average_ctd_profiles(
                yearMonth_o_c, indMali_o_c, tf_o_c)

        yearMonth_e_c, indMali_e_c, tf_e_c = _average_ctd_profiles(
                yearMonth_e_c, indMali_e_c, tf_e_c)

        print(f"tf_o_c: {tf_o_c.shape}")
        print(f"tf_e_c: {tf_e_c.shape}")

        if np.all(np.isnan(tf_o_c)):
            print(f"tf_o_c is all NaNs")
        
        if np.all(np.isnan(tf_e_c)):
            print(f"tf_e_c is all NaNs")

        # Distill MALI into corresponding yearMonth and cells for each orig and extrap
        ym_to_ind = {ym: i for i, ym in enumerate(self.yearMonth_mali)}
        ind_time_o = np.array([ym_to_ind[ym] for ym in yearMonth_o_c])
        tf_o_m = self.tf_mali[ind_time_o, indMali_o_c, :]
        print(f"tf_o_m: {tf_o_m.shape}")

        if np.all(np.isnan(tf_o_m)):
            print(f"tf_o_m is all NaNs")
        
        ind_time_e = np.array([ym_to_ind[ym] for ym in yearMonth_e_c])
        tf_e_m = self.tf_mali[ind_time_e, indMali_e_c, :]
        
        # remove profiles that are all zeros. 
        # <<TO DO>>: understand why this even happens
        #empty = np.all(tf_e_m == 0, axis=1)
        #tf_e_m = tf_e_m[~empty]
        #tf_e_c = tf_e_c[~empty]

        print(f"tf_e_m: {tf_e_m.shape}")

        if np.all(np.isnan(tf_e_m)):
            print(f"tf_e_m is all NaNs")
        
        # Compute skill score
        def _skill_score(model, obs):
            # remove all-NaN profiles
            valid = ~np.all(np.isnan(obs), axis=1)
            model = model[valid, :]
            obs = obs[valid, :]

            obs_mean = np.nanmean(obs, axis=1)
            num = np.nansum((model - obs) ** 2, axis=1)
            den = np.nansum((np.abs(model - obs_mean[:, None]) + np.abs(obs - obs_mean[:, None])) ** 2, axis=1)
            
            return 1.0 - num / den
            
        print("computing ECCO WSS")
        t0 = perf_counter()
        WSS_ecco = np.nanmean(_skill_score(tf_o_m, tf_o_c))
        print(f"WSS_ecco took {perf_counter() - t0:.1f} s")
        
        print("computing Extrap WSS")
        t0 = perf_counter()
        WSS_extrap = np.nanmean(_skill_score(tf_e_m, tf_e_c))
        print(f"WSS_extrap took {perf_counter() - t0:.1f} s")

        print(f"ECCO WSS: {WSS_ecco}")
        print(f"Extrap WSS: {WSS_extrap}")
        
        # save profiles
        ds = xr.Dataset()
        ds['tf_o_m'] = xr.DataArray(tf_o_m.astype('float64'), dims=("origProfiles","depth"))
        ds['tf_o_c'] = xr.DataArray(tf_o_c.astype('float64'), dims=("origProfiles","depth"))
        ds['tf_e_m'] = xr.DataArray(tf_e_m.astype('float64'), dims=("extrapProfiles","depth"))
        ds['tf_e_c'] = xr.DataArray(tf_e_c.astype('float64'), dims=("extrapProfiles","depth"))
        ds['z'] = xr.DataArray(self.z.astype('float64'), dims=("depth"))
        ds.to_netcdf('profiles.nc')


def _calc_thermal_forcing(temp, salt, depth):
    
    # from Jenkins 2011
    gamma1 = -5.7e-2
    gamma2 = 8.8e-2
    gamma3 = 7.61e-4

    freezingTemperature = (gamma1 * salt + gamma2 + gamma3 * depth)
    thermalForcing = temp - freezingTemperature

    return thermalForcing

def main():
    run = validateMaliExtrap()
    
    # Hadley EN4 is primary comparison dataset. Add others as desired
    run.process_Hadley_EN4()

    run.calc_willmott_skill_score()
 
if __name__ == "__main__":
    main()


