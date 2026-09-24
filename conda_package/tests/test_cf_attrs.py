import xarray

from mpas_tools.mesh.attrs import MESH_VAR_ATTRS, cf_conventions
from mpas_tools.mesh.conversion import convert, cull
from mpas_tools.planar_hex import make_planar_hex_mesh

from .util import get_test_data_file


def _check_mesh_attrs(ds):
    """
    Check that every mesh variable in the dataset has the attributes from
    ``MESH_VAR_ATTRS`` and that the file follows CF and MPAS conventions
    """
    for name in ds.data_vars:
        if name not in MESH_VAR_ATTRS:
            continue
        for key, value in MESH_VAR_ATTRS[name].items():
            assert ds[name].attrs.get(key) == value, (name, key)
        if 'units' not in MESH_VAR_ATTRS[name]:
            assert 'units' not in ds[name].attrs, name
    assert ds.attrs['Conventions'].split()[0].startswith('CF-')
    assert 'MPAS' in ds.attrs['Conventions'].split()


def test_cf_conventions():
    assert cf_conventions() == 'CF-1.8 MPAS'
    assert cf_conventions('') == 'CF-1.8 MPAS'
    assert cf_conventions('MPAS') == 'CF-1.8 MPAS'
    assert cf_conventions('CF-1.8 MPAS') == 'CF-1.8 MPAS'
    assert cf_conventions('MPAS CF-1.10') == 'CF-1.10 MPAS'
    assert cf_conventions('CF-1.6, ACDD-1.3') == 'CF-1.6 ACDD-1.3 MPAS'


def test_planar_hex_cf_attrs():
    ds = make_planar_hex_mesh(
        nx=10, ny=10, dc=1e3, nonperiodic_x=False, nonperiodic_y=True
    )
    _check_mesh_attrs(ds)
    assert ds.attrs['Conventions'] == 'CF-1.8 MPAS'


def test_cull_convert_planar_cf_attrs():
    ds = make_planar_hex_mesh(
        nx=10, ny=10, dc=1e3, nonperiodic_x=False, nonperiodic_y=True
    )
    ds_culled = cull(ds)
    _check_mesh_attrs(ds_culled)
    ds_converted = convert(ds_culled)
    _check_mesh_attrs(ds_converted)
    assert ds_converted.attrs['Conventions'] == 'CF-1.8 MPAS'


def test_convert_spherical_cf_attrs():
    ds = xarray.open_dataset(get_test_data_file('mesh.QU.1920km.151026.nc'))
    # an existing CF version should be kept
    ds.attrs['Conventions'] = 'MPAS CF-1.10'
    ds_converted = convert(ds)
    _check_mesh_attrs(ds_converted)
    assert ds_converted.attrs['Conventions'] == 'CF-1.10 MPAS'
    for name in ['cellQuality', 'gridSpacing', 'triangleQuality']:
        assert name in ds_converted
