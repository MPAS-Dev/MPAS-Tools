import os
import subprocess

import numpy as np
import pytest
import xarray as xr

from mpas_tools.io import write_netcdf

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
    ds = xr.open_dataset(TEST_MESH)
    out_file = tmp_path / 'test_cdf5.nc'
    write_netcdf(ds, str(out_file), format='NETCDF3_64BIT_DATA')
    # Use ncdump -k to check format
    result = subprocess.run(
        ['ncdump', '-k', str(out_file)],
        capture_output=True,
        text=True,
        check=True,
    )
    # Should be cdf5 for NETCDF3_64BIT_DATA
    assert result.stdout.strip() == 'cdf5'


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
