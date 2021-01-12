#!/usr/bin/env python
import numpy
import xarray
import os

from mpas_tools.io import write_netcdf
from mpas_tools.ocean.depth import compute_depth, compute_zmid, add_depth, \
    add_zmid, write_time_varying_zmid


def create_3d_mesh():
    outFileName = 'test_depth_mesh.nc'
    if os.path.exists(outFileName):
        dsMesh = xarray.open_dataset(outFileName)
    else:
        dsMesh = xarray.open_dataset(
            'mesh_tools/mesh_conversion_tools/test/mesh.QU.1920km.151026.nc')
        nCells = dsMesh.sizes['nCells']
        nVertLevels = 10
        zmax = 1000.
        layerThickness = zmax/nVertLevels
        dsMesh['refBottomDepth'] = \
            ('nVertLevels', numpy.linspace(layerThickness, zmax, nVertLevels))
        dsMesh['maxLevelCell'] = \
            ('nCells', nVertLevels*numpy.ones(nCells, dtype=int))
        dsMesh['bottomDepth'] = ('nCells', zmax*numpy.ones(nCells))
        dsMesh['layerThickness'] = \
            (('Time', 'nCells', 'nVertLevels'),
             layerThickness*numpy.ones((1, nCells, nVertLevels)))
        write_netcdf(dsMesh, 'test_depth_mesh.nc')

    return dsMesh


def test_compute_depth():
    dsMesh = create_3d_mesh()
    depth, depth_bnds = compute_depth(dsMesh.refBottomDepth)
    assert numpy.all(numpy.isclose(depth, numpy.linspace(50., 950., 10)))
    assert numpy.all(numpy.isclose(depth_bnds[:, 0],
                                   numpy.linspace(0., 900., 10)))
    assert numpy.all(numpy.isclose(depth_bnds[:, 1],
                                   numpy.linspace(100., 1000., 10)))


def test_compute_zmid():
    dsMesh = create_3d_mesh()
    zMid = compute_zmid(dsMesh.bottomDepth, dsMesh.maxLevelCell,
                        dsMesh.layerThickness, depth_dim='nVertLevels')

    assert zMid.dims == ('Time', 'nCells', 'nVertLevels')

    depth = zMid.isel(Time=0, nCells=0)
    assert numpy.all(numpy.isclose(depth.values,
                                   numpy.linspace(-50., -950., 10)))


def test_add_depth():
    dsMesh = create_3d_mesh()
    mesh_filename = 'test_depth_mesh.nc'
    out_filename = 'test_depth_out.nc'

    dsIn = xarray.Dataset()
    dsIn['temperature'] = xarray.ones_like(dsMesh.layerThickness)
    write_netcdf(dsIn, 'test_depth_in.nc')

    # test adding depth coordinate once to the mesh and once to the input file,
    # with the mesh passed in separately
    for in_filename, coord_filename in [(mesh_filename, None),
                                        ('test_depth_in.nc', mesh_filename)]:
        add_depth(in_filename, out_filename, coordFileName=coord_filename)
        dsOut = xarray.open_dataset(out_filename)
        assert 'depth' in dsOut.dims

        depth = dsOut.depth
        assert numpy.all(numpy.isclose(depth.values,
                                       numpy.linspace(50., 950., 10)))


def test_add_zmid():
    dsMesh = create_3d_mesh()
    mesh_filename = 'test_depth_mesh.nc'
    out_filename = 'test_depth_out.nc'

    dsIn = xarray.Dataset()
    dsIn['temperature'] = xarray.ones_like(dsMesh.layerThickness)
    write_netcdf(dsIn, 'test_depth_in.nc')

    # test adding zMid once to the mesh and once to the input file, with the
    # mesh passed in separately
    for in_filename, coord_filename in [(mesh_filename, None),
                                        ('test_depth_in.nc', mesh_filename)]:
        add_zmid(in_filename, out_filename, coordFileName=coord_filename)
        dsOut = xarray.open_dataset(out_filename)
        assert 'depth' in dsOut.dims

        zMid = dsOut.zMid
        assert zMid.dims == ('nCells', 'depth')

        depth = zMid.isel(nCells=0)
        assert numpy.all(numpy.isclose(depth.values,
                                       numpy.linspace(-50., -950., 10)))


def test_write_time_varying_zmid():

    dsMesh = create_3d_mesh()
    nCells = dsMesh.sizes['nCells']
    nVertLevels = dsMesh.sizes['nVertLevels']
    mesh_filename = 'test_depth_mesh.nc'
    in_filename = 'test_depth_in.nc'
    out_filename = 'test_depth_out.nc'

    layerThickness = 100.

    # test adding zMid once to the mesh and once to the input file, with the
    # mesh passed in separately, each one without and once with a prefix
    for coord_filename, prefix in [(None, ''),
                                   (mesh_filename, ''),
                                   (None, 'timeMonthly_avg_'),
                                   (mesh_filename, 'timeMonthly_avg_')]:
        print(coord_filename, prefix)

        if coord_filename is None:
            dsIn = dsMesh.drop_vars('layerThickness')
        else:
            dsIn = xarray.Dataset()
        layerThicknessVar = '{}layerThickness'.format(prefix)
        dsIn[layerThicknessVar] = \
            (('Time', 'nCells', 'nVertLevels'),
             layerThickness*numpy.ones((2, nCells, nVertLevels)))
        dsIn['{}temperature'.format(prefix)] = \
            xarray.ones_like(dsIn[layerThicknessVar])
        write_netcdf(dsIn, in_filename)
        dsIn.close()

        write_time_varying_zmid(in_filename, out_filename,
                                coordFileName=coord_filename, prefix=prefix)

        dsOut = xarray.open_dataset(out_filename)
        assert 'depth' in dsOut.dims

        zMid = dsOut['{}zMid'.format(prefix)]
        assert zMid.dims == ('Time', 'nCells', 'depth')

        depth = zMid.isel(Time=0, nCells=0)
        assert numpy.all(numpy.isclose(depth.values,
                                       numpy.linspace(-50., -950., 10)))
        dsOut.close()


if __name__ == '__main__':
    test_compute_depth()
    test_compute_zmid()
    test_add_depth()
    test_add_zmid()
    test_write_time_varying_zmid()
