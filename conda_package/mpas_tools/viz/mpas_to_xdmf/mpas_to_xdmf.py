from mpas_tools.viz.mpas_to_xdmf.io import _convert_to_xdmf, _load_dataset
from mpas_tools.viz.mpas_to_xdmf.mesh import _get_ds_mesh
from mpas_tools.viz.mpas_to_xdmf.time import _set_time


class MpasToXdmf:
    """
    A class for converting MPAS NetCDF files to XDMF + HDF5 format for
    visualization in ParaView.

    Attributes
    ----------
    ds : xarray.Dataset
        An xarray Dataset to be converted.
    ds_mesh : xarray.Dataset
        An xarray Dataset representing the mesh.
    """

    def __init__(self, ds=None, ds_mesh=None, xtime_var=None):
        """
        Initialize the converter with a mesh file and optional time series
        files.

        Parameters
        ----------
        ds : xarray.Dataset, optional
            An xarray Dataset, all variables from which will be converted. If
            ds_mesh is not provided, ds must also contain the mesh variables
            if it is provided.
        ds_mesh : xarray.Dataset, optional
            An xarray Dataset representing the mesh. These variables will not
            be converted (unless they are also in ds).
        xtime_var : str, optional
            Name of the variable containing time information, only used if
            ``ds`` is not ``None``.
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
        Load the MPAS mesh file and optionally combine it with time series
        files into a single xarray Dataset.

        Parameters
        ----------
        mesh_filename : str
            Path to the MPAS mesh file.
        time_series_filenames : list of str or str, optional
            List of filenames or a wildcard string for time series files.
        variables : list of str, optional
            List of variables to convert.
        xtime_var : str, optional
            Name of the variable containing time information.
        """
        self.ds_mesh, self.ds = _load_dataset(
            mesh_filename=mesh_filename,
            time_series_filenames=time_series_filenames,
            variables=variables,
            xtime_var=xtime_var,
        )

    def convert_to_xdmf(self, out_dir):
        """
        Convert an xarray Dataset to XDMF + HDF5 format.

        Parameters
        ----------
        out_dir : str
            Directory where XDMF and HDF5 files will be saved.
        """

        _convert_to_xdmf(
            ds=self.ds,
            ds_mesh=self.ds_mesh,
            out_dir=out_dir,
        )


def main():
    """
    Command-line interface for the MpasToXdmf.
    """
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

    args = parser.parse_args()

    converter = MpasToXdmf()
    converter.load(
        mesh_filename=args.mesh,
        time_series_filenames=args.time_series,
        variables=args.variables,
        xtime_var=args.xtime,
    )
    converter.convert_to_xdmf(out_dir=args.output_dir)
