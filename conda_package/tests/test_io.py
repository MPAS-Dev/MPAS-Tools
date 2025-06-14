import os
import subprocess

import numpy as np
import pytest
import xarray as xr

from mpas_tools.io import write_netcdf
from mpas_tools.logging import LoggingContext

from .util import get_test_data_file

TEST_MESH = get_test_data_file('mesh.QU.1920km.151026.nc')


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_write_netcdf_basic(tmp_path):
    ds = xr.open_dataset(TEST_MESH)
    out_file = tmp_path / 'test_basic.nc'
    write_netcdf(ds, str(out_file))
    ds2 = xr.open_dataset(out_file)
    # Should have same dimensions and variables
    assert set(ds.dims) == set(ds2.dims)
    for var in ds.data_vars:
        assert var in ds2.data_vars
    ds2.close()


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_write_netcdf_cdf5_format(tmp_path):
    with LoggingContext('test_write_netcdf_cdf5_format') as logger:
        ds = xr.open_dataset(TEST_MESH)
        out_file = tmp_path / 'test_cdf5.nc'
        # test with and without logging
        for logger_arg in [None, logger]:
            write_netcdf(
                ds,
                str(out_file),
                format='NETCDF3_64BIT_DATA',
                logger=logger_arg,
            )
    # Use ncdump -k to check format
    result = subprocess.run(
        ['ncdump', '-k', str(out_file)],
        capture_output=True,
        text=True,
        check=True,
    )
    # Should be cdf5 for NETCDF3_64BIT_DATA
    assert result.stdout.strip() == 'cdf5'
    # Check that the temporary file was deleted
    tmp_file = (
        out_file.parent / f'_tmp_{out_file.stem}.netcdf4{out_file.suffix}'
    )
    assert not os.path.exists(tmp_file)


def test_write_netcdf_int64_conversion_and_attr(tmp_path):
    # Create a dataset with int64 variable and an attribute
    arr = np.array([1, 2, 3], dtype=np.int64)
    ds = xr.Dataset({'foo': (('x',), arr)})
    ds['foo'].attrs['myattr'] = 'testattr'
    out_file = tmp_path / 'test_int64.nc'
    write_netcdf(ds, str(out_file))
    ds2 = xr.open_dataset(out_file)
    # Should be int32, not int64
    assert ds2['foo'].dtype == np.int32
    # Attribute should be preserved
    assert ds2['foo'].attrs['myattr'] == 'testattr'
    ds2.close()


def test_write_netcdf_fill_value(tmp_path):
    # Test that NaN values are written with correct fill value
    arr = np.array([1.0, np.nan, 3.0], dtype=np.float32)
    ds = xr.Dataset({'bar': (('x',), arr)})
    out_file = tmp_path / 'test_fill.nc'
    write_netcdf(ds, str(out_file))
    ds2 = xr.open_dataset(out_file)
    # The second value should be the default fill value for float32
    fill_value = ds2['bar'].encoding.get('_FillValue', None)
    assert fill_value is not None
    assert np.isnan(ds2['bar'].values[1])
    ds2.close()


def test_write_netcdf_string_dim_name(tmp_path):
    # Test that custom char_dim_name is used in encoding
    arr = np.array([b'abc', b'def'])
    ds = xr.Dataset({'baz': (('x',), arr)})
    out_file = tmp_path / 'test_strdim.nc'
    write_netcdf(ds, str(out_file), char_dim_name='CustomStrLen')
    ds2 = xr.open_dataset(out_file)
    # Should have the variable and correct shape
    assert 'baz' in ds2.variables
    ds2.close()
    arr = np.array([1, 2, 3], dtype=np.int64)
    ds = xr.Dataset({'foo': (('x',), arr)})
    ds['foo'].attrs['myattr'] = 'testattr'
    out_file = tmp_path / 'test_int64.nc'
    write_netcdf(ds, str(out_file))
    ds2 = xr.open_dataset(out_file)
    # Should be int32, not int64
    assert ds2['foo'].dtype == np.int32
    # Attribute should be preserved
    assert ds2['foo'].attrs['myattr'] == 'testattr'
    ds2.close()
