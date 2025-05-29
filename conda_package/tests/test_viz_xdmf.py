import os
import sys

import numpy as np
import pytest
import xarray as xr

from mpas_tools.io import write_netcdf
from mpas_tools.viz.mpas_to_xdmf.io import (
    _load_dataset,
    _parse_indices,
    _process_extra_dims,
)
from mpas_tools.viz.mpas_to_xdmf.mpas_to_xdmf import MpasToXdmf, main
from mpas_tools.viz.mpas_to_xdmf.time import _set_time

from .util import get_test_data_file

TEST_MESH = get_test_data_file('mesh.QU.1920km.151026.nc')


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_load_mesh_only():
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)
    assert isinstance(converter.ds, xr.Dataset)
    assert isinstance(converter.ds_mesh, xr.Dataset)
    # Should have mesh dimensions
    assert 'nCells' in converter.ds.dims


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_set_time_with_no_xtime():
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)
    # Should create a 'Time' variable if 'Time' in dims
    if 'Time' in converter.ds.dims:
        assert 'Time' in converter.ds
        arr = converter.ds['Time'].values
        assert np.all(arr == np.arange(converter.ds.sizes['Time']))


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_convert_to_xdmf(tmp_path):
    converter = MpasToXdmf()
    variables = ['xCell', 'areaCell', 'cellsOnCell']
    extra_dims = {'maxEdges': [0]}
    converter.load(mesh_filename=TEST_MESH, variables=variables)
    out_dir = tmp_path / 'out'
    converter.convert_to_xdmf(str(out_dir), extra_dims=extra_dims)
    # Check that output files exist for cells
    assert (out_dir / 'fieldsOnCells.h5').exists()
    assert (out_dir / 'fieldsOnCells.xdmf').exists()


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_extra_dims(tmp_path):
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)
    # Simulate an extra dimension if present
    extra_dims = {}
    for dim in converter.ds.dims:
        if dim not in ['Time', 'nCells', 'nEdges', 'nVertices']:
            extra_dims[dim] = [0]
    out_dir = tmp_path / 'out_extra'
    converter.convert_to_xdmf(str(out_dir), extra_dims=extra_dims)
    assert (out_dir / 'fieldsOnCells.h5').exists()


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_load_with_time_series_and_variables(tmp_path):
    ts1 = tmp_path / 'ts1.nc'
    ts2 = tmp_path / 'ts2.nc'

    # Simulate a time series by adding xtime and area variables
    ds = xr.open_dataset(TEST_MESH)
    ds['xtime'] = ('Time', ['0001-01-01_00:00:00'])
    ds['area'] = (('Time', 'nCells'), ds.areaCell.values[None, :])
    write_netcdf(ds, ts1)
    ds['xtime'] = ('Time', ['0001-01-02_00:00:00'])
    write_netcdf(ds, ts2)

    variables = ['areaCell', 'area']

    converter = MpasToXdmf()
    converter.load(
        mesh_filename=TEST_MESH,
        time_series_filenames=[str(ts1), str(ts2)],
        variables=variables,
    )
    print(converter.ds)
    for var in variables:
        assert var in converter.ds.data_vars, (
            f'Variable {var} not found in dataset'
        )
    assert converter.ds.sizes['Time'] == 2


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_process_extra_dims_drop(tmp_path):
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)

    # drop all variables with extra dimensions
    extra_dims = {
        'maxEdges': [],
        'maxEdges2': [],
        'TWO': [],
        'vertexDegree': [],
    }

    ds = _process_extra_dims(converter.ds, extra_dims=extra_dims)
    for dim in extra_dims:
        assert dim not in ds.dims, f'Dimension {dim} should be dropped'


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_set_time_invalid_xtime(tmp_path):
    ts1 = tmp_path / 'ts1.nc'
    # Simulate a time-depndent variable and add xtime
    ds = xr.open_dataset(TEST_MESH)
    ds['xtime'] = ('Time', ['0001-01-01_00:00:00'])
    ds['area'] = (('Time', 'nCells'), ds.areaCell.values[None, :])
    write_netcdf(ds, ts1)

    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH, time_series_filenames=[str(ts1)])
    # Should raise ValueError if xtime_var is not present
    with pytest.raises(ValueError):
        _set_time(ds=converter.ds, xtime_var='not_a_var')


def test_parse_indices_invalid_cases():
    # Should raise on mixed slice/list
    with pytest.raises(ValueError):
        _parse_indices('1:3,5', 5)
    # Should raise on invalid string
    with pytest.raises(ValueError):
        _parse_indices('foo', 5)


def test_parse_indices_valid_cases():
    # Empty list
    assert _parse_indices('', 5) == []
    # Single index
    assert _parse_indices('0', 5) == [0]
    # Comma-separated list
    assert _parse_indices('1,2,3', 5) == [1, 2, 3]
    # Slice notation
    assert _parse_indices('0:3', 5) == [0, 1, 2]
    # Slice with stride
    assert _parse_indices('0:5:2', 5) == [0, 2, 4]
    # Full slice
    assert _parse_indices(':', 4) == [0, 1, 2, 3]


def test_main_cli(monkeypatch, tmp_path):
    # Test CLI entry point with minimal arguments
    mesh = TEST_MESH
    if not os.path.exists(mesh):
        pytest.skip('Test mesh not available')
    out_dir = tmp_path / 'cli_out'
    sys_argv = ['prog', '-m', mesh, '-o', str(out_dir), '-v', 'areaCell']
    monkeypatch.setattr(sys, 'argv', sys_argv)
    # Patch input to always return blank (skip extra dims)
    monkeypatch.setattr('builtins.input', lambda _: '')
    main()
    assert (out_dir / 'fieldsOnCells.h5').exists()


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_load_dataset_missing_variable():
    # Should not raise if variable is missing in mesh, but should raise if not
    # present at all
    with pytest.raises(KeyError):
        _load_dataset(TEST_MESH, None, ['not_a_var'], None)
