from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import netCDF4
from datetime import datetime
import sys


def write_netcdf(ds, fileName, fillValues=netCDF4.default_fillvals,
                 format='NETCDF3_64BIT'):
    '''Write an xarray Dataset with NetCDF4 fill values where needed'''
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
