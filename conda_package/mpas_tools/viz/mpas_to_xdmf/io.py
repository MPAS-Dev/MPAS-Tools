import importlib.resources
import os

import h5py
import xarray as xr
from jinja2 import Template

from mpas_tools.viz.mpas_to_xdmf.geometry import (
    _build_cell_geometry,
    _build_edge_geometry,
    _build_vertex_geometry,
)
from mpas_tools.viz.mpas_to_xdmf.mesh import _get_ds_mesh
from mpas_tools.viz.mpas_to_xdmf.time import _set_time


def _load_dataset(mesh_filename, time_series_filenames, variables, xtime_var):
    """
    Load the MPAS mesh file and optionally combine it with time series
    files into a mesh dataset and a dataset to convert.

    Parameters
    ----------
    mesh_filename : str
        Path to the MPAS mesh file.
    time_series_filenames : list of str or str
        List of filenames or a wildcard string for time series files.
    variables : list of str
        List of variables to convert. Special keys include:
        - 'allOnCells': Load all variables with dimension `nCells`.
        - 'allOnEdges': Load all variables with dimension `nEdges`.
        - 'allOnVertices': Load all variables with dimension `nVertices`.
    xtime_var : str
        Name of the variable containing time information.

    Returns
    -------
    ds_mesh : xarray.Dataset
        The mesh dataset.
    ds : xarray.Dataset
        The dataset to convert, which may include the mesh variables.
    """
    # Load the mesh file
    ds_mesh = xr.open_dataset(mesh_filename)

    if time_series_filenames is None:
        ds = ds_mesh
    else:
        ds = xr.open_mfdataset(
            time_series_filenames,
            combine='nested',
            concat_dim='Time',
            data_vars='minimal',
            coords='minimal',
            compat='override',
        )

    ds_mesh = _get_ds_mesh(ds_mesh)
    if variables is not None:
        selected_vars = set()
        for var in variables:
            if var == 'allOnCells':
                selected_vars.update(
                    [v for v in ds.data_vars if 'nCells' in ds[v].dims]
                )
            elif var == 'allOnEdges':
                selected_vars.update(
                    [v for v in ds.data_vars if 'nEdges' in ds[v].dims]
                )
            elif var == 'allOnVertices':
                selected_vars.update(
                    [v for v in ds.data_vars if 'nVertices' in ds[v].dims]
                )
            else:
                selected_vars.add(var)
        ds = ds[list(selected_vars)]

    _set_time(ds=ds, xtime_var=xtime_var)

    return ds_mesh, ds


def _convert_to_xdmf(ds, ds_mesh, out_dir):
    """
    Convert an xarray Dataset to XDMF + HDF5 format.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to convert.
    ds_mesh : xarray.Dataset
        The mesh dataset.
    out_dir : str
        Directory where XDMF and HDF5 files will be saved.
    """
    os.makedirs(out_dir, exist_ok=True)

    if 'nCells' in ds.dims:
        _convert_cells_to_xdmf(ds, ds_mesh, out_dir)
    if 'nEdges' in ds.dims:
        _convert_edges_to_xdmf(ds, ds_mesh, out_dir)
    if 'nVertices' in ds.dims:
        _convert_vertices_to_xdmf(ds, ds_mesh, out_dir)


def _convert_cells_to_xdmf(ds, ds_mesh, out_dir):
    """
    Convert cell-centered data to XDMF + HDF5 format.
    """
    ds_cell_geom = _build_cell_geometry(ds_mesh)
    cell_vars = [var for var in ds.data_vars if 'nCells' in ds[var].dims]
    ds_cells = ds[cell_vars]
    _write_xdmf(ds_cell_geom, ds_cells, out_dir, suffix='Cells')


def _convert_edges_to_xdmf(ds, ds_mesh, out_dir):
    """
    Convert edge-centered data to XDMF + HDF5 format.
    """
    ds_edge_geom = _build_edge_geometry(ds_mesh)
    edge_vars = [var for var in ds.data_vars if 'nEdges' in ds[var].dims]
    ds_edges = ds[edge_vars]
    _write_xdmf(ds_edge_geom, ds_edges, out_dir, suffix='Edges')


def _convert_vertices_to_xdmf(ds, ds_mesh, out_dir):
    """
    Convert vertex-centered data to XDMF + HDF5 format.
    """
    ds_vertex_geom = _build_vertex_geometry(ds_mesh)
    vertex_vars = [var for var in ds.data_vars if 'nVertices' in ds[var].dims]
    ds_vertices = ds[vertex_vars]
    _write_xdmf(ds_vertex_geom, ds_vertices, out_dir, suffix='Vertices')


def _write_xdmf(ds_geom, ds_data, out_dir, suffix):
    """
    Write data to HDF5 and metadata to XDMF format.

    Parameters
    ----------
    ds_geom : xarray.Dataset
        Dataset containing geometry information (e.g., points, connectivity).
    ds_data : xarray.Dataset
        Dataset containing time-varying data to write.
    out_dir : str
        Directory where XDMF and HDF5 files will be saved.
    suffix : str
        Suffix to append to output filenames (e.g., 'Cells', 'Edges').
    """
    h5_basename = f'fieldsOn{suffix}.h5'
    h5_filename = os.path.join(out_dir, h5_basename)
    xdmf_filename = os.path.join(out_dir, f'fieldsOn{suffix}.xdmf')

    # Write HDF5 file
    with h5py.File(h5_filename, 'w') as h5_file:
        # Write geometry
        h5_file.create_dataset('Points', data=ds_geom['points'].values)
        h5_file.create_dataset('Cells', data=ds_geom['cells'].values)

        # Write time-varying and static data
        for var_name in ds_data.data_vars:
            if 'Time' in ds_data[var_name].dims:
                for t_idx in range(ds_data.sizes['Time']):
                    dataset_name = f'{var_name}_t{t_idx}'
                    da = ds_data[var_name].isel(Time=t_idx)
                    h5_file.create_dataset(dataset_name, data=da.values)
            else:
                h5_file.create_dataset(var_name, data=ds_data[var_name].values)

    # Preprocess variable metadata for the template
    variables_metadata = [
        {
            'name': var_name,
            'has_time': 'Time' in ds_data[var_name].dims,
        }
        for var_name in ds_data.data_vars
    ]

    # Load XDMF template from external file
    package = 'mpas_tools.viz.mpas_to_xdmf.templates'
    filename = 'xdmf_template.xml'
    with importlib.resources.open_text(package, filename) as template_file:
        xdmf_template = Template(template_file.read())

    # Render XDMF file
    cells = ds_geom['cells'].values
    times = ds_data['Time'].values if 'Time' in ds_data.dims else []

    xdmf_content = xdmf_template.render(
        times=times,
        num_elements=cells.shape[0],
        num_verts=cells.shape[1],
        num_points=ds_geom['points'].shape[0],
        variables=variables_metadata,
        suffix=suffix,
        h5_basename=h5_basename,
    )

    with open(xdmf_filename, 'w') as xdmf_file:
        xdmf_file.write(xdmf_content)
