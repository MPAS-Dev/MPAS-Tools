"""
Extract Cartesian (X, Y, Z), zonal and meridional components of an MPAS
vector field, given the field on edge normals.

This tool requires that the field 'coeffs_reconstruct' has been saved to a
NetCDF file.  The simplest way to do this is to include the following
stream in a forward run:

<stream name="vector_reconstruction"
        clobber_mode="truncate"
        type="output"
        output_interval="initial_only"
        filename_template="vector_reconstruction.nc">

    <var name="coeffs_reconstruct"/>
</stream>

and run the model for one time step.
"""
import argparse
import sys
from datetime import datetime

import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

from mpas_tools.io import write_netcdf


def reconstruct_variable(out_var_name, variable_on_edges, ds_mesh,
                         coeffs_reconstruct, ds_out, chunk_size=32768,
                         quiet=False):
    """
    Extract Cartesian (X, Y, Z), zonal and meridional components of an MPAS
    vector field, given the field on edge normals.

    Parameters
    ----------
    out_var_name : str
        The output variable name

    variable_on_edges : xarray.DataArray
        The variable at edge normals

    ds_mesh : xarray.Dataset
        A dataset with the mesh variables

    coeffs_reconstruct : xarray.DataArray
        A data array with the reconstruction coefficients

    ds_out : xarray.Dataset
        The dataset the output variables should be added to

    chunk_size : int, optional
        The size of chunks along the ``nCells`` dimension used to divide up
        the work

    quiet : bool, optional
        Whether to suppress print statements and progress bars
    """
    n_cells = ds_mesh.sizes['nCells']
    edges_on_cell = ds_mesh.edgesOnCell - 1

    variable_on_edges.load()
    edges_on_cell.load()
    coeffs_reconstruct.load()

    dims = []
    sizes = []
    varIndices = {}
    for dim in variable_on_edges.dims:
        size = variable_on_edges.sizes[dim]
        varIndices[dim] = np.arange(size)
        if dim == 'nEdges':
            dim = 'nCells'
            size = n_cells
            varIndices['nEdges'] = edges_on_cell
        dims.append(dim)
        sizes.append(size)

    coeffs_reconstruct = coeffs_reconstruct.chunk({'nCells': chunk_size})

    variable = variable_on_edges[varIndices].chunk({'nCells': chunk_size})
    if quiet:
        variable.compute()
    else:
        print(f'Computing {out_var_name} at edgesOnCell:')
        with ProgressBar():
            variable.compute()

    var_cart = []

    if not quiet:
        print('Computing Cartesian components:')
    for index, component in enumerate(['X', 'Y', 'Z']):
        var = (coeffs_reconstruct.isel(R3=index)*variable).sum(
            dim='maxEdges').transpose(*dims)
        out_name = f'{out_var_name}{component}'
        if quiet:
            var.compute()
        else:
            print(out_name)
            with ProgressBar():
                var.compute()
        ds_out[out_name] = var
        var_cart.append(var)

    lat_cell = ds_mesh.latCell
    lon_cell = ds_mesh.lonCell
    lat_cell.load()
    lon_cell.load()

    clat = np.cos(lat_cell)
    slat = np.sin(lat_cell)
    clon = np.cos(lon_cell)
    slon = np.sin(lon_cell)

    if not quiet:
        print('Computing zonal and meridional components:')

    out_name = f'{out_var_name}Zonal'
    zonal = -var_cart[0] * slon + var_cart[1] * clon
    if quiet:
        zonal.compute()
    else:
        print(out_name)
        with ProgressBar():
            zonal.compute()
    ds_out[out_name] = zonal

    out_name = f'{out_var_name}Meridional'
    merid = (-(var_cart[0] * clon + var_cart[1] * slon) * slat +
             var_cart[2] * clat)
    if quiet:
        merid.compute()
    else:
        print(out_name)
        with ProgressBar():
            merid.compute()
    ds_out[out_name] = merid


def main():
    # client = Client(n_workers=1, threads_per_worker=4, memory_limit='10GB')
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-m", "--mesh_filename", dest="mesh_filename",
                        type=str, required=False,
                        help="An MPAS file with mesh data (edgesOnCell, etc.) "
                             "if not from in_filename")
    parser.add_argument("-w", "--weights_filename", dest="weights_filename",
                        type=str, required=False,
                        help="An MPAS file with coeffs_reconstruct if not "
                             "from in_filename")
    parser.add_argument("-i", "--in_filename", dest="in_filename", type=str,
                        required=True,
                        help="An MPAS file with one or more fields on edges "
                             "to be reconstructed at cell centers.  Used for "
                             "mesh data and/or weights if a separate files "
                             "are not provided.")
    parser.add_argument("-v", "--variables", nargs='+', dest="variables",
                        type=str, required=True,
                        help="variables on edges to reconstruct")
    parser.add_argument("--out_variables", nargs='+', dest="out_variables",
                        type=str, required=False,
                        help="prefixes for output variable names")
    parser.add_argument("-o", "--out_filename", dest="out_filename", type=str,
                        required=True,
                        help="An output MPAS file with the reconstructed "
                             "X, Y, Z, zonal and meridional fields")
    args = parser.parse_args()

    if args.mesh_filename:
        mesh_filename = args.mesh_filename
    else:
        mesh_filename = args.in_filename

    if args.weights_filename:
        weights_filename = args.weights_filename
    else:
        weights_filename = args.in_filename

    if args.out_variables is not None:
        out_variables = args.out_variables
    else:
        out_variables = args.variables

    ds_in = xr.open_dataset(args.in_filename, mask_and_scale=False)
    ds_mesh = xr.open_dataset(mesh_filename)
    ds_weights = xr.open_dataset(weights_filename)
    coeffs_reconstruct = ds_weights.coeffs_reconstruct
    ds_out = xr.Dataset()

    for in_var_name, out_var_name in zip(args.variables, out_variables):
        reconstruct_variable(out_var_name, ds_in[in_var_name], ds_mesh,
                             coeffs_reconstruct, ds_out)

    for attr_name in ds_in.attrs:
        ds_out.attrs[attr_name] = ds_in.attrs[attr_name]

    time = datetime.now().strftime('%c')

    history = f'{time}: {" ".join(sys.argv)}'

    if 'history' in ds_out.attrs:
        ds_out.attrs['history'] = f'{history}\n{ds_out.attrs["history"]}'
    else:
        ds_out.attrs['history'] = history

    write_netcdf(ds_out, args.out_filename)


if __name__ == '__main__':
    main()
