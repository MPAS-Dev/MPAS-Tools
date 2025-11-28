"""
MPAS to XDMF Converter
======================

This module provides the :py:class:`MpasToXdmf` class and a command-line
interface for converting MPAS NetCDF files to XDMF + HDF5 format, suitable
for visualization in ParaView and other tools.

Features
--------
- Supports cell-, edge-, and vertex-centered variables.
- Allows slicing along extra dimensions (e.g., nVertLevels).
- Combines mesh and time series files.
- Selects variables by name or by special keys (e.g., 'allOnCells').

Example Usage
-------------
Python:

    >>> from mpas_tools.viz.mpas_to_xdmf.mpas_to_xdmf import MpasToXdmf
    >>> converter = MpasToXdmf()
    >>> converter.load(mesh_filename="mesh.nc", time_series_filenames="output.*.nc",
    ...                variables=["temperature", "salinity"], xtime_var="xtime")
    >>> converter.convert_to_xdmf(out_dir="output_dir", extra_dims={"nVertLevels": [0, 1, 2]})

Command line:

    $ mpas_to_xdmf -m mesh.nc -t output.*.nc -v temperature salinity -o output_dir -d nVertLevels=0:3

See Also
--------
- Documentation: https://mpas-dev.github.io/MPAS-Tools/latest/mpas_to_xdmf.html
"""  # noqa: E501

from mpas_tools.viz.mpas_to_xdmf.io import (
    _convert_to_xdmf,
    _load_dataset,
    _parse_extra_dims,
    _process_extra_dims,
)
from mpas_tools.viz.mpas_to_xdmf.mesh import _get_ds_mesh
from mpas_tools.viz.mpas_to_xdmf.time import _set_time


class MpasToXdmf:
    """
    A class for converting MPAS NetCDF files to XDMF + HDF5 format for
    visualization in ParaView.

    Attributes
    ----------
    ds : xarray.Dataset
        The dataset containing variables to convert.
    ds_mesh : xarray.Dataset
        The mesh dataset (may be the same as ds if no time series is provided).

    Notes
    -----
    - Use the `load()` method to read mesh and time series files.
    - Use `convert_to_xdmf()` to write XDMF and HDF5 files.
    - Special variable keys: 'allOnCells', 'allOnEdges', 'allOnVertices'.
    - Extra dimensions can be sliced using the `extra_dims` argument.
    """

    def __init__(self, ds=None, ds_mesh=None, xtime_var=None):
        """
        Initialize the converter with a mesh file and optional time series
        files.

        Parameters
        ----------
        ds : xarray.Dataset, optional
            An xarray Dataset containing variables to convert. If ds_mesh is
            not provided, ds must also contain mesh variables.
        ds_mesh : xarray.Dataset, optional
            An xarray Dataset representing the mesh. If not provided, ds is
            used as the mesh.
        xtime_var : str, optional
            Name of the variable containing time information (e.g., 'xtime').
        """
        if ds is not None and ds_mesh is None:
            ds_mesh = ds
        if ds_mesh is not None:
            self.ds_mesh = _get_ds_mesh(ds_mesh)
        self.ds = ds
        if ds is not None:
            _set_time(ds=ds, xtime_var=xtime_var)

    def load(
        self,
        mesh_filename,
        time_series_filenames=None,
        variables=None,
        xtime_var=None,
    ):
        """
        Load the MPAS mesh file and optional time series.

        Parameters
        ----------
        mesh_filename : str
            Path to the MPAS mesh file (NetCDF).
        time_series_filenames : list of str or str, optional
            List of NetCDF filenames or a wildcard string for time series
            files. If None, only the mesh file is used.
        variables : list of str, optional
            List of variables to convert. Special keys:

            * ``"allOnCells"``: all variables with dimension ``"nCells"``.
            * ``"allOnEdges"``: all variables with dimension ``"nEdges"``.
            * ``"allOnVertices"``: all variables with dimension
              ``"nVertices"``.

            If None, all variables are included.
        xtime_var : str, optional
            Name of the variable containing time information (e.g.,
            ``"xtime"``).
        """
        self.ds_mesh, self.ds = _load_dataset(
            mesh_filename=mesh_filename,
            time_series_filenames=time_series_filenames,
            variables=variables,
            xtime_var=xtime_var,
        )

    def convert_to_xdmf(self, out_dir, extra_dims=None, quiet=False):
        """
        Convert the loaded xarray Dataset to XDMF + HDF5 format.

        Parameters
        ----------
        out_dir : str
            Directory where XDMF and HDF5 files will be saved.
        extra_dims : dict, optional
            Dictionary mapping extra dimensions to their selected indices.
            Example - ``{'nVertLevels': [0, 1, 2]}``. If None, all indices are
            included.
        quiet : bool, optional
            If True, suppress progress output.

        Output
        ------
        - XDMF files (.xdmf) and HDF5 files (.h5) in the specified directory.

        Raises
        ------
        ValueError
            If required files or variables are missing.
        """
        # Process extra dimensions
        self.ds = _process_extra_dims(self.ds, extra_dims)

        _convert_to_xdmf(
            ds=self.ds,
            ds_mesh=self.ds_mesh,
            out_dir=out_dir,
            quiet=quiet,
        )


def main():
    """
    Command-line interface for the MpasToXdmf.

    Usage
    -----
    $ mpas_to_xdmf -m mesh.nc -t output.*.nc -v temperature salinity -o output_dir -d nVertLevels=0:3

    See `mpas_to_xdmf --help` for all options.
    """  # noqa: E501
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            'Convert MPAS NetCDF files to XDMF + HDF5 format for ParaView.'
        )
    )
    parser.add_argument(
        '-m', '--mesh', required=True, help='Path to the MPAS mesh file.'
    )
    parser.add_argument(
        '-t',
        '--time-series',
        nargs='*',
        help=(
            'Wildcard or list of time series files.  If not provided, the '
            'mesh file will be used.'
        ),
    )
    parser.add_argument(
        '-v',
        '--variables',
        nargs='*',
        help=(
            'List of variables to convert. Special keys include: '
            "'allOnCells' (variables with dimension `nCells`), "
            "'allOnEdges' (variables with dimension `nEdges`), and "
            "'allOnVertices' (variables with dimension `nVertices`). "
            'If not provided, all variables will be converted.'
        ),
    )
    parser.add_argument(
        '-o',
        '--output-dir',
        required=True,
        help='Directory to save XDMF files.',
    )
    parser.add_argument(
        '-x',
        '--xtime',
        help='Name of the variable containing time information.',
    )
    parser.add_argument(
        '-d',
        '--dim-list',
        nargs='+',
        help=(
            'List of dimensions and indices to slice '
            '(e.g., nVertLevels=0:10:2).'
        ),
    )
    parser.add_argument(
        '-q',
        '--quiet',
        action='store_true',
        help='Suppress progress output.',
    )

    args = parser.parse_args()

    converter = MpasToXdmf()
    converter.load(
        mesh_filename=args.mesh,
        time_series_filenames=args.time_series,
        variables=args.variables,
        xtime_var=args.xtime,
    )

    # Parse extra dimensions
    extra_dims = _parse_extra_dims(args.dim_list, converter.ds)

    converter.convert_to_xdmf(
        out_dir=args.output_dir,
        extra_dims=extra_dims,
        quiet=args.quiet,
    )
