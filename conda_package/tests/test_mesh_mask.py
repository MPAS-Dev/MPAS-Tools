#!/usr/bin/env python

import multiprocessing

import numpy as np
import pyproj
import xarray as xr
from geometric_features import FeatureCollection, GeometricFeatures

from mpas_tools.cime.constants import constants
from mpas_tools.io import write_netcdf
from mpas_tools.mesh.conversion import convert, cull
from mpas_tools.mesh.mask import (
    compute_lon_lat_region_masks,
    compute_mpas_flood_fill_mask,
    compute_mpas_region_masks,
    compute_mpas_transect_masks,
    compute_projection_grid_region_masks,
)


def test_compute_mpas_region_masks():
    ds_mesh, _ = _get_mesh()
    gf = GeometricFeatures()
    fc_mask = gf.read(
        componentName='ocean',
        objectType='region',
        featureNames=['North Atlantic Ocean'],
    )

    pool = _get_pool()
    ds_masks = compute_mpas_region_masks(
        ds_mesh, fc_mask, maskTypes=('cell', 'edge', 'vertex'), pool=pool
    )

    for dim in ['nCells', 'nEdges', 'nVertices', 'nRegions']:
        assert dim in ds_masks.dims

    for var in [
        'regionCellMasks',
        'regionEdgeMasks',
        'regionVertexMasks',
        'regionNames',
    ]:
        assert var in ds_masks.data_vars


def test_compute_mpas_transect_masks_no_edge_sign():
    ds_mesh, earth_radius = _get_mesh()
    gf = GeometricFeatures()
    fc_mask = gf.read(
        componentName='ocean',
        objectType='transect',
        featureNames=['Drake Passage'],
    )

    pool = _get_pool()
    ds_masks = compute_mpas_transect_masks(
        ds_mesh,
        fc_mask,
        earth_radius,
        maskTypes=('cell', 'edge', 'vertex'),
        pool=pool,
        addEdgeSign=False,
    )
    # write_netcdf(ds_masks, 'transect_masks.nc')

    for dim in ['nCells', 'nEdges', 'nVertices', 'nTransects']:
        assert dim in ds_masks.dims

    for var in [
        'transectCellMasks',
        'transectEdgeMasks',
        'transectVertexMasks',
        'transectNames',
    ]:
        assert var in ds_masks.data_vars


def test_compute_mpas_transect_masks_edge_sign():
    ds_mesh, earth_radius = _get_mesh()
    # write_netcdf(ds_mesh, 'mesh.nc')
    gf = GeometricFeatures()
    fc_mask = gf.read(
        componentName='ocean',
        objectType='transect',
        featureNames=['Drake Passage'],
    )

    pool = _get_pool()
    ds_masks = compute_mpas_transect_masks(
        ds_mesh,
        fc_mask,
        earth_radius,
        maskTypes=('edge',),
        pool=pool,
        addEdgeSign=True,
    )
    # write_netcdf(ds_masks, 'transect_edge_mask.nc')

    for dim in ['nEdges', 'nTransects']:
        assert dim in ds_masks.dims

    for var in ['transectEdgeMasks', 'transectEdgeMaskSigns', 'transectNames']:
        assert var in ds_masks.data_vars


def test_compute_mpas_flood_fill_mask():
    ds_mesh, earth_radius = _get_mesh()
    gf = GeometricFeatures()
    fc_mask = gf.read(
        componentName='ocean',
        objectType='region',
        featureNames=['Global Ocean 15S to 15N'],
    )

    pool = _get_pool()
    ds_mask = compute_mpas_region_masks(
        ds_mesh, fc_mask, maskTypes=('cell',), pool=pool
    )

    ds_mesh = cull(ds_mesh, ds_mask)
    # write_netcdf(ds_mesh, 'culled_mesh.nc')

    feature = {
        'type': 'Feature',
        'properties': {
            'name': 'Mid North Atlantic',
            'tags': '',
            'object': 'point',
            'component': 'ocean',
            'author': 'Xylar Asay-Davis',
        },
        'geometry': {'type': 'Point', 'coordinates': [-40.000000, 40.000000]},
    }
    fc_seed = FeatureCollection()
    fc_seed.add_feature(feature)

    ds_mask = compute_mpas_flood_fill_mask(ds_mesh, fc_seed)
    # write_netcdf(ds_mask, 'mask_floodfill.nc')

    assert 'nCells' in ds_mask.dims
    assert 'cellSeedMask' in ds_mask.data_vars


def test_compute_lon_lat_region_masks():
    lon = np.linspace(-180, 180, 37)
    lat = np.linspace(-90, 90, 19)
    gf = GeometricFeatures()
    fc_mask = gf.read(
        componentName='ocean',
        objectType='region',
        featureNames=['North Atlantic Ocean'],
    )
    pool = _get_pool()
    ds_mask = compute_lon_lat_region_masks(lon, lat, fc_mask, pool=pool)
    # write_netcdf(ds_mask, 'mask_lon_lat.nc')

    for dim in ['lon', 'lat', 'nRegions']:
        assert dim in ds_mask.dims

    for var in ['regionMasks', 'regionNames']:
        assert var in ds_mask.data_vars


def test_compute_projection_grid_region_masks():
    x = np.linspace(-6000e3, 6000e3, 121)
    y = np.linspace(-6000e3, 6000e3, 121)

    projection = pyproj.Proj(
        '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
        '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'
    )
    lat_lon_projection = pyproj.Proj(proj='latlong', datum='WGS84')

    transformer = pyproj.Transformer.from_proj(projection, lat_lon_projection)
    x_grid, y_grid = np.meshgrid(x, y)
    lon, lat = transformer.transform(x_grid, y_grid)
    gf = GeometricFeatures()
    fc_mask = gf.read(
        componentName='landice',
        objectType='region',
        featureNames=['ISMIP6 Basin A-Ap'],
    )
    pool = _get_pool()
    ds_mask = compute_projection_grid_region_masks(
        lon, lat, fc_mask, pool=pool
    )
    ds_mask['lon'] = (('y', 'x'), lon)
    ds_mask['lat'] = (('y', 'x'), lat)
    write_netcdf(ds_mask, 'mask_proj.nc')

    for dim in ['x', 'y', 'nRegions']:
        assert dim in ds_mask.dims

    for var in ['regionMasks', 'regionNames']:
        assert var in ds_mask.data_vars


def _get_pool():
    multiprocessing.set_start_method('forkserver', force=True)
    process_count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(process_count)
    return pool


def _get_mesh():
    ds_mesh = xr.open_dataset(
        'mesh_tools/mesh_conversion_tools/test/mesh.QU.1920km.151026.nc'
    )
    earth_radius = constants['SHR_CONST_REARTH']
    ds_mesh.attrs['sphere_radius'] = earth_radius
    for coord in [
        'xCell',
        'yCell',
        'zCell',
        'xVertex',
        'yVertex',
        'zVertex',
        'xEdge',
        'yEdge',
        'zEdge',
        'dcEdge',
        'dvEdge',
    ]:
        ds_mesh[coord] = earth_radius * ds_mesh[coord]
    ds_mesh = convert(ds_mesh)
    return ds_mesh, earth_radius


if __name__ == '__main__':
    test_compute_mpas_region_masks()
    test_compute_mpas_transect_masks_no_edge_sign()
    test_compute_mpas_transect_masks_edge_sign()
    test_compute_mpas_flood_fill_mask()
    test_compute_lon_lat_region_masks()
    test_compute_projection_grid_region_masks()
