import glob
import importlib.resources
import os

import h5py
import xarray as xr
from jinja2 import Template
from tqdm import tqdm

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
        if isinstance(time_series_filenames, str):
            # Deterministic ordering for reproducibility
            file_list = sorted(glob.glob(time_series_filenames))
        else:
            file_list = time_series_filenames

        if len(file_list) == 0:
            raise ValueError(
                f'No time series files found matching '
                f"'{time_series_filenames}'"
            )

        if len(file_list) == 1:
            # If only one file, open it directly
            ds = xr.open_dataset(file_list[0])
        else:
            ds = xr.open_mfdataset(
                file_list,
                combine='nested',
                concat_dim='Time',
                data_vars='minimal',
                coords='minimal',
                compat='override',
                decode_times=False,
                decode_timedelta=False,
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
        if xtime_var is not None:
            selected_vars.add(xtime_var)
        ds = ds[list(selected_vars)]

    _set_time(ds=ds, xtime_var=xtime_var)

    return ds_mesh, ds


def _convert_to_xdmf(ds, ds_mesh, out_dir, quiet=False):
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
    quiet : bool, optional
        If True, suppress progress output. Default is False.
    """
    os.makedirs(out_dir, exist_ok=True)

    if 'nCells' in ds.dims:
        _convert_cells_to_xdmf(ds, ds_mesh, out_dir, quiet)
    if 'nEdges' in ds.dims:
        _convert_edges_to_xdmf(ds, ds_mesh, out_dir, quiet)
    if 'nVertices' in ds.dims:
        _convert_vertices_to_xdmf(ds, ds_mesh, out_dir, quiet)


def _convert_cells_to_xdmf(ds, ds_mesh, out_dir, quiet):
    """
    Convert cell-centered data to XDMF + HDF5 format.
    """
    ds_cell_geom = _build_cell_geometry(ds_mesh)
    cell_vars = [var for var in ds.data_vars if 'nCells' in ds[var].dims]
    ds_cells = ds[cell_vars]
    _write_xdmf(ds_cell_geom, ds_cells, out_dir, suffix='Cells', quiet=quiet)


def _convert_edges_to_xdmf(ds, ds_mesh, out_dir, quiet):
    """
    Convert edge-centered data to XDMF + HDF5 format.
    """
    ds_edge_geom = _build_edge_geometry(ds_mesh)
    edge_vars = [var for var in ds.data_vars if 'nEdges' in ds[var].dims]
    ds_edges = ds[edge_vars]
    _write_xdmf(ds_edge_geom, ds_edges, out_dir, suffix='Edges', quiet=quiet)


def _convert_vertices_to_xdmf(ds, ds_mesh, out_dir, quiet):
    """
    Convert vertex-centered data to XDMF + HDF5 format.
    """
    ds_vertex_geom = _build_vertex_geometry(ds_mesh)
    vertex_vars = [var for var in ds.data_vars if 'nVertices' in ds[var].dims]
    vert_to_kite_map = ds_vertex_geom['vert_to_kite_map']
    ds_vertices = ds[vertex_vars].isel(nVertices=vert_to_kite_map)
    _write_xdmf(
        ds_vertex_geom, ds_vertices, out_dir, suffix='Vertices', quiet=quiet
    )


def _write_xdmf(ds_geom, ds_data, out_dir, suffix, quiet=False):
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
    quiet : bool, optional
        If True, suppress progress output. Default is False.
    """
    h5_basename = f'fieldsOn{suffix}.h5'
    h5_filename = os.path.join(out_dir, h5_basename)
    xdmf_filename = os.path.join(out_dir, f'fieldsOn{suffix}.xdmf')

    # Write HDF5 file
    with h5py.File(h5_filename, 'w') as h5_file:
        # Write geometry
        h5_file.create_dataset('Points', data=ds_geom['points'].values)
        h5_file.create_dataset('Cells', data=ds_geom['cells'].values)

        # Calculate total progress steps
        total_steps = sum(
            ds_data.sizes['Time'] if 'Time' in ds_data[var].dims else 1
            for var in ds_data.data_vars
        )

        # Write time-varying and static data with progress bar
        if quiet:
            iterator = None
        else:
            iterator = tqdm(
                ds_data.data_vars, total=total_steps, desc=f'Writing {suffix}'
            )
        for var_name in ds_data.data_vars:
            if iterator is not None:
                iterator.set_description(f'Processing {var_name}')
            if 'Time' in ds_data[var_name].dims:
                for t_idx in range(ds_data.sizes['Time']):
                    dataset_name = f'{var_name}_t{t_idx}'
                    da = ds_data[var_name].isel(Time=t_idx)
                    h5_file.create_dataset(dataset_name, data=da.values)
                    if iterator is not None:
                        iterator.update(1)
            else:
                h5_file.create_dataset(var_name, data=ds_data[var_name].values)
                if iterator is not None:
                    iterator.update(1)

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


def _parse_extra_dims(dimension_list, ds):
    """
    Parse and prompt for indices of extra dimensions.

    Parameters
    ----------
    dimension_list : list of str
        List of dimensions and indices in the format <dimension>=<indices>.
    ds : xarray.Dataset
        Dataset to extract dimensions from.

    Returns
    -------
    extra_dims : dict
        Dictionary mapping dimensions to their selected indices.
    """
    extra_dims = {}
    unspecified_dims = []

    for dim in ds.dims:
        if dim in ['Time', 'nCells', 'nEdges', 'nVertices']:
            continue

        # Check if the dimension is specified in the dimension_list
        specified = False
        if dimension_list:
            for dim_spec in dimension_list:
                if dim_spec.startswith(f'{dim}='):
                    indices = dim_spec.split('=')[1]
                    extra_dims[dim] = _parse_indices(indices, ds.sizes[dim])
                    specified = True
                    break

        if not specified:
            unspecified_dims.append(dim)

    # If there are unspecified dimensions, display a detailed prompt
    if unspecified_dims:
        print('\nThe following dimensions require indices to be specified:')
        print('You can enter indices in one of the following formats:')
        print("  - A single index (e.g., '0')")
        print("  - A comma-separated list of indices (e.g., '0,1,2')")
        print("  - A range with optional stride (e.g., ':' or '0:10:2')")
        print('  - Leave blank to skip fields with this dimension.')
        print()

        for dim in unspecified_dims:
            print(f"Dimension '{dim}' has size {ds.sizes[dim]}.")
            indices = input(f"Enter indices to keep for '{dim}': ")
            extra_dims[dim] = _parse_indices(indices, ds.sizes[dim])

    return extra_dims


def _parse_indices(index_string, dim_size):
    """
    Parse an index string into a list of indices.

    Parameters
    ----------
    index_string : str
        Index string (e.g., "0,1,2", "0:10:2", ":").
    dim_size : int
        Size of the dimension.

    Returns
    -------
    indices : list of int
        Parsed indices.
    """
    if not index_string:
        return []
    if ':' in index_string:
        # Support slice notation like ':', '0:10', '0:10:2', etc.
        parts = index_string.split(':')
        # Validate that parts has at most 3 elements
        if len(parts) > 3:
            raise ValueError(
                f"Invalid index string '{index_string}': too many colons. "
                'Expected at most two colons.'
            )
        # Pad parts to length 3 with empty strings if needed
        while len(parts) < 3:
            parts.append('')
        # Convert to int or None
        start = int(parts[0]) if parts[0] else 0
        stop = int(parts[1]) if parts[1] else dim_size
        step = int(parts[2]) if parts[2] else 1
        return list(range(start, stop, step))
    return [int(i) for i in index_string.split(',')]


def _process_extra_dims(ds, extra_dims):
    """
    Process extra dimensions in the dataset by ensuring all are covered,
    unwrapping variables with extra dimensions into multiple variables with
    basic dimensions, and applying slicing or dropping variables as needed.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to process for extra dimensions.
    extra_dims : dict
        Dictionary mapping dimensions to their selected indices.

    Returns
    -------
    ds : xarray.Dataset
        The processed dataset with extra dimensions handled.

    Raises
    ------
    ValueError
        If any extra dimensions are not covered by the extra_dims dictionary.
    """
    basic_dims = ['Time', 'nCells', 'nEdges', 'nVertices']
    for dim in ds.dims:
        if dim not in basic_dims and dim not in (extra_dims or {}):
            raise ValueError(
                f"Dimension '{dim}' is not covered by the extra_dims "
                f'dictionary.'
            )

    if extra_dims:
        for dim, indices in extra_dims.items():
            if not indices:
                # Drop variables with the given dimension if the list is empty
                ds = ds.drop_vars(
                    [var for var in ds.data_vars if dim in ds[var].dims]
                )
            else:
                # Unwrap variables with the extra dimension
                vars_to_unwrap = [
                    var for var in ds.data_vars if dim in ds[var].dims
                ]
                for var in vars_to_unwrap:
                    for index in indices:
                        # Create a new variable with the suffix `_index`
                        new_var_name = f'{var}_{index}'
                        ds[new_var_name] = ds[var].isel({dim: index})
                    # Drop the original variable with the extra dimension
                    ds = ds.drop_vars(var)

    return ds
