#!/usr/bin/env python
import xarray as xr

from mpas_tools.cime.constants import constants
from mpas_tools.mesh.conversion import convert, cull
from mpas_tools.planar_hex import make_planar_hex_mesh
from mpas_tools.translate import center
from mpas_tools.viz.mesh_to_triangles import mesh_to_triangles
from mpas_tools.viz.transect import (
    find_planar_transect_cells_and_weights,
    find_spherical_transect_cells_and_weights,
    make_triangle_tree,
)


def test_mesh_to_triangles():
    _, ds_tris = _get_triangles()
    for dim in ['nTriangles', 'nNodes', 'nInterp']:
        assert dim in ds_tris.dims

    for var in [
        'triCellIndices',
        'nodeCellIndices',
        'nodeCellWeights',
        'xNode',
        'yNode',
        'zNode',
        'latNode',
        'lonNode',
    ]:
        assert var in ds_tris.data_vars


def test_make_triangle_tree():
    _, ds_tris = _get_triangles()
    tree = make_triangle_tree(ds_tris)
    assert tree.n == ds_tris.sizes['nTriangles'] * ds_tris.sizes['nNodes']


def test_find_spherical_transect_cells_and_weights():
    ds_mesh, ds_tris = _get_triangles()
    tree = make_triangle_tree(ds_tris)
    lon_transect = xr.DataArray([-10.0, 0.0, 10.0], dims=('nPoints',))
    lat_transect = xr.DataArray([-10.0, 0.0, 10.0], dims=('nPoints',))
    ds_transect = find_spherical_transect_cells_and_weights(
        lon_transect,
        lat_transect,
        ds_tris,
        ds_mesh,
        tree,
        degrees=True,
        subdivision_res=10e3,
    )

    for dim in [
        'nSegments',
        'nNodes',
        'nHorizBounds',
        'nHorizWeights',
        'nPoints',
    ]:
        assert dim in ds_transect.dims

    for var in [
        'xCartNode',
        'yCartNode',
        'zCartNode',
        'dNode',
        'lonNode',
        'latNode',
        'horizTriangleIndices',
        'horizCellIndices',
        'horizTriangleNodeIndices',
        'interpHorizCellIndices',
        'interpHorizCellWeights',
        'lonTransect',
        'latTransect',
        'xCartTransect',
        'yCartTransect',
        'zCartTransect',
        'dTransect',
        'transectIndicesOnHorizNode',
        'transectWeightsOnHorizNode',
    ]:
        assert var in ds_transect.data_vars


def test_find_planar_transect_cells_and_weights():
    ds_mesh = make_planar_hex_mesh(
        nx=102, ny=52, dc=4e3, nonperiodic_x=True, nonperiodic_y=True
    )
    ds_mesh = cull(ds_mesh)
    center(ds_mesh)
    ds_tris = mesh_to_triangles(ds_mesh)
    tree = make_triangle_tree(ds_tris)

    x_transect = 1e3 * xr.DataArray([-10.0, 0.0, 10.0], dims=('nPoints',))
    y_transect = 1e3 * xr.DataArray([-10.0, 0.0, 10.0], dims=('nPoints',))

    ds_transect = find_planar_transect_cells_and_weights(
        x_transect, y_transect, ds_tris, ds_mesh, tree, subdivision_res=1e3
    )

    for dim in [
        'nSegments',
        'nNodes',
        'nHorizBounds',
        'nHorizWeights',
        'nPoints',
    ]:
        assert dim in ds_transect.dims

    for var in [
        'xNode',
        'yNode',
        'dNode',
        'horizTriangleIndices',
        'horizCellIndices',
        'horizTriangleNodeIndices',
        'interpHorizCellIndices',
        'interpHorizCellWeights',
        'xTransect',
        'yTransect',
        'dTransect',
        'transectIndicesOnHorizNode',
        'transectWeightsOnHorizNode',
    ]:
        assert var in ds_transect.data_vars


def _get_triangles():
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
    ds_tris = mesh_to_triangles(ds_mesh)
    return ds_mesh, ds_tris


if __name__ == '__main__':
    test_mesh_to_triangles()
    test_make_triangle_tree()
    test_find_spherical_transect_cells_and_weights()
    test_find_planar_transect_cells_and_weights()
