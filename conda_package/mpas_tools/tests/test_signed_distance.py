#!/usr/bin/env python

import numpy as np
import xarray as xr

from geometric_features import FeatureCollection

from mpas_tools.cime.constants import constants
from mpas_tools.io import write_netcdf
from mpas_tools.mesh.creation.signed_distance import \
    signed_distance_from_geojson, distance_from_geojson, mask_from_geojson


def test_signed_distance_from_geojson():

    lon, lat, fc_mask, earth_radius = _get_lon_lat_fc()

    signed_distance = signed_distance_from_geojson(
        fc_mask, lon, lat, earth_radius, max_length=5.)

    ds = xr.Dataset()
    ds['lon'] = ('lon', lon)
    ds['lat'] = ('lat', lat)
    ds['signed_distance'] = (('lat', 'lon'), signed_distance)
    write_netcdf(ds, 'signed_distance.nc')


def test_distance_from_geojson():
    lon, lat, fc_mask, earth_radius = _get_lon_lat_fc()

    distance = distance_from_geojson(
        fc_mask, lon, lat, earth_radius, max_length=5.)

    ds = xr.Dataset()
    ds['lon'] = ('lon', lon)
    ds['lat'] = ('lat', lat)
    ds['distance'] = (('lat', 'lon'), distance)
    write_netcdf(ds, 'distance.nc')


def test_mask_from_geojson():
    lon, lat, fc_mask, earth_radius = _get_lon_lat_fc()

    mask = mask_from_geojson(fc_mask, lon, lat)

    ds = xr.Dataset()
    ds['lon'] = ('lon', lon)
    ds['lat'] = ('lat', lat)
    ds['mask'] = (('lat', 'lon'), mask)
    write_netcdf(ds, 'mask.nc')


def _get_lon_lat_fc():

    lon = np.linspace(-180, 180, 37)
    lat = np.linspace(-90, 90, 19)
    earth_radius = constants['SHR_CONST_REARTH']

    feature = {
        "type": "Feature",
        "properties": {
            "name": "North Atlantic",
            "tags": "",
            "object": "region",
            "component": "ocean",
            "author": "Xylar Asay-Davis"
        },
        "geometry": {
            "type": "Polygon",
            "coordinates": [
                [
                    [
                        -39.53161633291441,
                        57.08649995213068
                    ],
                    [
                        -69.30597933223675,
                        28.03212363054105
                    ],
                    [
                        -33.914428822900334,
                        13.4287331666755
                    ],
                    [
                        -15.9735991802836,
                        39.731395665957876
                    ],
                    [
                        -39.53161633291441,
                        57.08649995213068
                    ]
                ]
            ]
        }
    }
    fc_mask = FeatureCollection()
    fc_mask.add_feature(feature)
    return lon, lat, fc_mask, earth_radius


if __name__ == '__main__':
    test_signed_distance_from_geojson()
    test_distance_from_geojson()
    test_mask_from_geojson()
