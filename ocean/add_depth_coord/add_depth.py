#!/usr/bin/env python

"""
Add a 1D coordinate "depth" to an MPAS-Ocean output file that defines the
positive-up vertical location of each layer.
"""
# Authors
# -------
# Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import numpy
import netCDF4
import argparse
import sys


def write_netcdf(ds, fileName, fillValues=netCDF4.default_fillvals):
    '''
    Write an xarray data set to a NetCDF file using finite fill values

    Parameters
    ----------
    ds : xarray.Dataset object
        The xarray data set to be written to a file

    fileName : str
        The fileName to write the data set to

    fillValues : dict
        A dictionary of fill values for each supported data type.  By default,
        this is the dictionary used by the netCDF4 package.  Key entries should
        be of the form 'f8' (for float64), 'i4' (for int32), etc.
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    encodingDict = {}
    variableNames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variableName in variableNames:
        dtype = ds[variableName].dtype
        for fillType in fillValues:
            if dtype == numpy.dtype(fillType):
                encodingDict[variableName] = \
                    {'_FillValue': fillValues[fillType]}
                break

    ds.to_netcdf(fileName, encoding=encodingDict)


def compute_depth(refBottomDepth):
    """
    Computes depth given refBottomDepth

    Parameters
    ----------
    refBottomDepth : ``xarray.DataArray``
        the depth of the bottom of each vertical layer in the initial state
        (perfect z-level coordinate)

    Returns
    -------
    depth : ``xarray.DataArray``
        the vertical coordinate defining the middle of each layer
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    refBottomDepth = refBottomDepth.values

    depth = numpy.zeros(refBottomDepth.shape)

    depth[0] = 0.5*refBottomDepth[0]
    depth[1:] = 0.5*(refBottomDepth[1:] + refBottomDepth[0:-1])

    return depth


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--coordFileName", dest="coordFileName",
                        type=str, required=False,
                        help="A MPAS-Ocean file with refBottomDepth")
    parser.add_argument("-i", "--inFileName", dest="inFileName", type=str,
                        required=True,
                        help="An input MPAS-Ocean file that depth should be"
                             "added to, used for coords if another file is"
                             "not provided via -c.")
    parser.add_argument("-o", "--outFileName", dest="outFileName", type=str,
                        required=True,
                        help="An output MPAS-Ocean file with depth added")
    args = parser.parse_args()

    if args.coordFileName:
        coordFileName = args.coordFileName
    else:
        coordFileName = args.inputFileName

    dsCoord = xarray.open_dataset(coordFileName)
    dsCoord = dsCoord.rename({'nVertLevels': 'depth'})

    ds = xarray.open_dataset(args.inFileName)
    ds = ds.rename({'nVertLevels': 'depth'})

    ds.coords['depth'] = ('depth',
                          compute_depth(dsCoord.refBottomDepth))
    ds.depth.attrs['units'] = 'meters'
    ds.depth.attrs['positive'] = 'down'
    ds.depth.attrs['standard_name'] = 'depth'

    ds.depth.attrs['long_name'] = 'reference depth of the center of each ' \
                                  'vertical level'

    for varName in ds.data_vars:
        var = ds[varName]
        if 'depth' in var.dims:
            var = var.assign_coords(depth=ds.depth)
            ds[varName] = var

    if 'history' in ds.attrs:
        ds.attrs['history'] = '{}\n{}'.format(' '.join(sys.argv),
                                              ds.attrs['history'])
    else:
        ds.attrs['history'] = ' '.join(sys.argv)

    write_netcdf(ds, args.outFileName)


if __name__ == '__main__':
    main()
