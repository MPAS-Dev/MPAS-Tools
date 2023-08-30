from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import netCDF4
from datetime import datetime
import sys


default_format = 'NETCDF3_64BIT'
default_engine = None
default_char_dim_name = 'StrLen'
default_fills = netCDF4.default_fillvals


def write_netcdf(ds, fileName, fillValues=None, format=None, engine=None,
                 char_dim_name=None):
    """
    Write an xarray.Dataset to a file with NetCDF4 fill values and the given
    name of the string dimension.  Also adds the time and command-line to the
    history attribute.

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

    """
    if format is None:
        format = default_format

    if fillValues is None:
        fillValues = default_fills

    if engine is None:
        engine = default_engine

    if char_dim_name is None:
        char_dim_name = default_char_dim_name

    encodingDict = {}
    variableNames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variableName in variableNames:
        isNumeric = numpy.issubdtype(ds[variableName].dtype, numpy.number)
        if isNumeric and numpy.any(numpy.isnan(ds[variableName])):
            dtype = ds[variableName].dtype
            for fillType in fillValues:
                if dtype == numpy.dtype(fillType):
                    encodingDict[variableName] = \
                        {'_FillValue': fillValues[fillType]}
                    break
        else:
            encodingDict[variableName] = {'_FillValue': None}

        isString = numpy.issubdtype(ds[variableName].dtype, numpy.string_)
        if isString and char_dim_name is not None:
            encodingDict[variableName] = {'char_dim_name': char_dim_name}

    update_history(ds)

    if 'Time' in ds.dims:
        # make sure the Time dimension is unlimited because MPAS has trouble
        # reading Time otherwise
        ds.encoding['unlimited_dims'] = {'Time'}

    ds.to_netcdf(fileName, encoding=encodingDict, format=format, engine=engine)


def update_history(ds):
    '''Add or append history to attributes of a data set'''

    thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + \
        " ".join(sys.argv[:])
    if 'history' in ds.attrs:
        newhist = '\n'.join([thiscommand, ds.attrs['history']])
    else:
        newhist = thiscommand
    ds.attrs['history'] = newhist
