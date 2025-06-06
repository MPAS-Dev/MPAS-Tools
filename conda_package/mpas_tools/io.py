import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import netCDF4
import numpy

from mpas_tools.logging import check_call

default_format = 'NETCDF3_64BIT'
default_engine = None
default_char_dim_name = 'StrLen'
default_fills = netCDF4.default_fillvals


def write_netcdf(
    ds,
    fileName,
    fillValues=None,
    format=None,
    engine=None,
    char_dim_name=None,
    logger=None,
):
    """
    Write an xarray.Dataset to a file with NetCDF4 fill values and the given
    name of the string dimension.  Also adds the time and command-line to the
    history attribute.

    Note: the ``NETCDF3_64BIT_DATA`` format is handled as a special case
    because xarray output with this format is not performant. First, the file
    is written in `NETCDF4` format, which supports larger files and variables.
    Then, the `ncks` command is used to convert the file to the
    `NETCDF3_64BIT_DATA` format.

    Note: All int64 variables are automatically converted to int32 for MPAS
    compatibility.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to save

    fileName : str
        The path for the NetCDF file to write

    fillValues : dict, optional
        A dictionary of fill values for different NetCDF types.  Default is
        ``mpas_tools.io.default_fills``, which can be modified but which
        defaults to ``netCDF4.default_fillvals``

    format : {'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC'}, optional
        The NetCDF file format to use.  Default is
        ``mpas_tools.io.default_format``, which can be modified but which
        defaults to ``'NETCDF3_64BIT'``

    engine : {'netcdf4', 'scipy', 'h5netcdf'}, optional
        The library to use for NetCDF output.  The default is the same as
        in :py:meth:`xarray.Dataset.to_netcdf` and depends on ``format``.
        You can override the default by setting
        ``mpas_tools.io.default_engine``

    char_dim_name : str, optional
        The name of the dimension used for character strings, or None to let
        xarray figure this out. Default is
        ``mpas_tools.io.default_char_dim_name``, which can be modified but
        which defaults to ``'StrLen'``

    logger : logging.Logger, optional
        A logger to write messages to write the output of `ncks` conversion
        calls to.  If None, `ncks` output is suppressed.  This is only
        relevant if `format` is 'NETCDF3_64BIT_DATA'
    """  # noqa: E501
    if format is None:
        format = default_format

    if fillValues is None:
        fillValues = default_fills

    if engine is None:
        engine = default_engine

    if char_dim_name is None:
        char_dim_name = default_char_dim_name

    # Convert int64 variables to int32 for MPAS compatibility
    for var in list(ds.data_vars.keys()) + list(ds.coords.keys()):
        if ds[var].dtype == numpy.int64:
            attrs = ds[var].attrs.copy()
            ds[var] = ds[var].astype(numpy.int32)
            ds[var].attrs = attrs

    encodingDict = {}
    variableNames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variableName in variableNames:
        isNumeric = numpy.issubdtype(ds[variableName].dtype, numpy.number)
        if isNumeric and numpy.any(numpy.isnan(ds[variableName])):
            dtype = ds[variableName].dtype
            for fillType in fillValues:
                if dtype == numpy.dtype(fillType):
                    encodingDict[variableName] = {
                        '_FillValue': fillValues[fillType]
                    }
                    break
        else:
            encodingDict[variableName] = {'_FillValue': None}

        isString = numpy.issubdtype(ds[variableName].dtype, numpy.bytes_)
        if isString and char_dim_name is not None:
            encodingDict[variableName] = {'char_dim_name': char_dim_name}

    update_history(ds)

    if 'Time' in ds.dims:
        # make sure the Time dimension is unlimited because MPAS has trouble
        # reading Time otherwise
        ds.encoding['unlimited_dims'] = {'Time'}

    # for performance, we have to handle this as a special case
    convert = format == 'NETCDF3_64BIT_DATA'

    if convert:
        out_path = Path(fileName)
        out_filename = (
            out_path.parent / f'_tmp_{out_path.stem}.netcdf4{out_path.suffix}'
        )
        format = 'NETCDF4'
        if engine == 'scipy':
            # that's not going to work
            engine = 'netcdf4'
    else:
        out_filename = fileName

    ds.to_netcdf(
        out_filename, encoding=encodingDict, format=format, engine=engine
    )

    if convert:
        args = [
            'ncks',
            '-O',
            '-5',
            out_filename,
            fileName,
        ]
        if logger is None:
            subprocess.run(
                args,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        else:
            check_call(args, logger=logger)
        # delete the temporary NETCDF4 file
        os.remove(out_filename)


def update_history(ds):
    """Add or append history to attributes of a data set"""

    thiscommand = (
        datetime.now().strftime('%a %b %d %H:%M:%S %Y')
        + ': '
        + ' '.join(sys.argv[:])
    )
    if 'history' in ds.attrs:
        newhist = '\n'.join([thiscommand, ds.attrs['history']])
    else:
        newhist = thiscommand
    ds.attrs['history'] = newhist
