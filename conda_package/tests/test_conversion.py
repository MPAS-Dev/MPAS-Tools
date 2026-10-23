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


def test_conversion_triangle_angle_quality():
    ds_mesh = xarray.open_dataset(
        get_test_data_file('mesh.QU.1920km.151026.nc')
    )
    ds_mesh = convert(dsIn=ds_mesh)

    edges_on_vertex = ds_mesh.edgesOnVertex.values - 1
    assert np.all(edges_on_vertex >= 0)
    dc_edge = ds_mesh.dcEdge.values
    a_len = dc_edge[edges_on_vertex[:, 0]]
    b_len = dc_edge[edges_on_vertex[:, 1]]
    c_len = dc_edge[edges_on_vertex[:, 2]]

    # law of cosines for the angle opposite each side of the dual triangle
    angle1 = np.arccos(
        np.clip((b_len**2 + c_len**2 - a_len**2) / (2 * b_len * c_len), -1, 1)
    )
    angle2 = np.arccos(
        np.clip((a_len**2 + c_len**2 - b_len**2) / (2 * a_len * c_len), -1, 1)
    )
    angle3 = np.arccos(
        np.clip((a_len**2 + b_len**2 - c_len**2) / (2 * a_len * b_len), -1, 1)
    )
    angles = np.stack([angle1, angle2, angle3], axis=1)
    assert np.max(np.abs(angles.sum(axis=1) - np.pi)) < 1.0e-10

    min_angle = angles.min(axis=1)
    max_angle = angles.max(axis=1)
    quality = ds_mesh.triangleAngleQuality.values
    assert np.max(np.abs(quality - min_angle / max_angle)) < 1.0e-10

    obtuse = (max_angle > 0.5 * np.pi).astype(int)
    assert np.array_equal(ds_mesh.obtuseTriangle.values, obtuse)


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
    test_conversion_triangle_angle_quality()
