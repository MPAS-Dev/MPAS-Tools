import xarray

from mpas_tools.mesh.attrs import (
    SPHERE_VAR_ATTRS,
    add_mesh_attrs,
    cf_conventions,
    get_mesh_var_attrs,
)
from mpas_tools.mesh.conversion import convert, cull
from mpas_tools.planar_hex import make_planar_hex_mesh

from .util import get_test_data_file


def _check_mesh_attrs(ds):
    """
    Check that every mesh variable in the dataset has the attributes from
    ``get_mesh_var_attrs()`` and that the file follows CF and MPAS conventions
    """
    on_a_sphere = ds.attrs['on_a_sphere'].strip() == 'YES'
    mesh_var_attrs = get_mesh_var_attrs(on_a_sphere)
    for name in ds.data_vars:
        if name not in mesh_var_attrs:
            continue
        for key, value in mesh_var_attrs[name].items():
            assert ds[name].attrs.get(key) == value, (name, key)
        if 'units' not in mesh_var_attrs[name]:
            assert 'units' not in ds[name].attrs, name
        if 'standard_name' not in mesh_var_attrs[name]:
            assert 'standard_name' not in ds[name].attrs, name
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


def test_get_mesh_var_attrs():
    for name, attrs in SPHERE_VAR_ATTRS.items():
        assert 'standard_name' not in get_mesh_var_attrs(False)[name]
        assert (
            get_mesh_var_attrs(True)[name]['standard_name']
            == attrs['standard_name']
        )


def test_add_mesh_attrs_planar_lat_lon():
    ds = make_planar_hex_mesh(
        nx=4, ny=4, dc=10e3, nonperiodic_x=True, nonperiodic_y=True
    )
    for name in SPHERE_VAR_ATTRS:
        assert 'standard_name' not in ds[name].attrs, name
    # planar meshes made with older versions have these standard names
    for name, attrs in SPHERE_VAR_ATTRS.items():
        ds[name].attrs.update(attrs)
    add_mesh_attrs(ds)
    _check_mesh_attrs(ds)


def test_add_mesh_attrs_spherical_lat_lon():
    ds = xarray.open_dataset(get_test_data_file('mesh.QU.1920km.151026.nc'))
    for name in SPHERE_VAR_ATTRS:
        ds[name].attrs = {}
    add_mesh_attrs(ds)
    _check_mesh_attrs(ds)
    assert ds.latCell.attrs['standard_name'] == 'latitude'
    assert ds.lonCell.attrs['standard_name'] == 'longitude'


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
