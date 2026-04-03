#!/usr/bin/env python

import matplotlib
import numpy as np

from mpas_tools.io import write_netcdf
from mpas_tools.mesh.conversion import convert, cull, mask
from mpas_tools.mesh.spherical import recompute_angle_edge

from .util import get_test_data_file

matplotlib.use('Agg')
import xarray
from geometric_features import read_feature_collection


def test_conversion():
    dsMesh = xarray.open_dataset(
        get_test_data_file('mesh.QU.1920km.151026.nc')
    )
    dsMesh = convert(dsIn=dsMesh)
    write_netcdf(dsMesh, 'mesh.nc')

    dsMask = xarray.open_dataset(get_test_data_file('land_mask_final.nc'))
    dsCulled = cull(dsIn=dsMesh, dsMask=dsMask)
    write_netcdf(dsCulled, 'culled_mesh.nc')

    fcMask = read_feature_collection(
        get_test_data_file('Arctic_Ocean.geojson')
    )
    dsMask = mask(dsMesh=dsMesh, fcMask=fcMask)
    write_netcdf(dsMask, 'antarctic_mask.nc')


def test_conversion_angle_edge():
    ds_mesh = xarray.open_dataset(
        get_test_data_file('mesh.QU.1920km.151026.nc')
    )
    ds_mesh = convert(dsIn=ds_mesh)

    angle_edge_python = recompute_angle_edge(ds_mesh)
    angle_diff = np.angle(
        np.exp(1j * (angle_edge_python.values - ds_mesh.angleEdge.values))
    )

    assert np.all(np.isfinite(angle_diff))
    assert np.max(np.abs(angle_diff)) < 1.0e-10


if __name__ == '__main__':
    test_conversion()
    test_conversion_angle_edge()
