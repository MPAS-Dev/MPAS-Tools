#!/usr/bin/env python

import matplotlib
import numpy as np

from mpas_tools.io import write_netcdf
from mpas_tools.mesh.conversion import _masks_to_int, convert, cull, mask
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


def test_masks_to_int_dataset_copy():
    ds_in = xarray.Dataset(
        data_vars={
            'regionCellMasks': (
                ('nCells', 'nRegions'),
                np.array([[True, False], [False, True]]),
            ),
            'cullCell': (('nCells',), np.array([False, True])),
            'xCell': (('nCells',), np.array([1.0, 2.0])),
        },
        attrs={'meshName': 'unit-test'},
    )

    ds_out = _masks_to_int(ds_in)

    assert ds_out.regionCellMasks.dtype == np.int32
    assert ds_out.cullCell.dtype == np.int32
    assert ds_out.attrs == ds_in.attrs
    assert np.array_equal(ds_out.xCell.values, ds_in.xCell.values)


if __name__ == '__main__':
    test_conversion()
    test_conversion_angle_edge()
