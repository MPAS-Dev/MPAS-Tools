#!/usr/bin/env python

from mpas_tools.mesh.conversion import convert, cull, mask
from mpas_tools.io import write_netcdf
import matplotlib
matplotlib.use('Agg')
from geometric_features import read_feature_collection
import xarray


def test_conversion():
    dsMesh = xarray.open_dataset(
        'mesh_tools/mesh_conversion_tools/test/mesh.QU.1920km.151026.nc')
    dsMesh = convert(dsIn=dsMesh)
    write_netcdf(dsMesh, 'mesh.nc')

    dsMask = xarray.open_dataset(
        'mesh_tools/mesh_conversion_tools/test/land_mask_final.nc')
    dsCulled = cull(dsIn=dsMesh, dsMask=dsMask)
    write_netcdf(dsCulled, 'culled_mesh.nc')

    dsMask = xarray.open_dataset(
        'mesh_tools/mesh_conversion_tools/test/land_mask_final.nc')

    fcMask = read_feature_collection(
        'mesh_tools/mesh_conversion_tools/test/Arctic_Ocean.geojson')
    dsMask = mask(dsMesh=dsMesh, fcMask=fcMask)
    write_netcdf(dsMask, 'antarctic_mask.nc')


if __name__ == '__main__':
    test_conversion()
