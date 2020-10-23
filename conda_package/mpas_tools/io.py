from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import netCDF4
from datetime import datetime
import sys


def write_netcdf(ds, fileName, fillValues=netCDF4.default_fillvals,
                 format='NETCDF3_64BIT', char_dim_name='StrLen'):
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
        A dictionary of fill values for different NetCDF types

    format : {'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT',
              'NETCDF3_CLASSIC'}, optional
        The NetCDF file format to use

    char_dim_name : str, optional
        The name of the dimension used for character strings, or None to let
        xarray figure this out.

    """
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

    ds.to_netcdf(fileName, encoding=encodingDict, format=format)


def update_history(ds):
    '''Add or append history to attributes of a data set'''

    thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + \
        " ".join(sys.argv[:])
    if 'history' in ds.attrs:
        newhist = '\n'.join([thiscommand, ds.attrs['history']])
    else:
        newhist = thiscommand
    ds.attrs['history'] = newhist
