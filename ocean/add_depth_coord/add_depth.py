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
import argparse
import sys
from datetime import datetime


def write_netcdf(ds, fileName):
    '''
    Write an xarray data set to a NetCDF file making use of the _FillValue
    attributes of each variable.  This function should be used for data sets
    opened with mask_and_scale=False.

    Parameters
    ----------
    ds : xarray.Dataset object
        The xarray data set to be written to a file

    fileName : str
        The fileName to write the data set to
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    encodingDict = {}
    variableNames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variableName in variableNames:
        if '_FillValue' in ds[variableName].attrs:
            encodingDict[variableName] = \
                {'_FillValue': ds[variableName].attrs['_FillValue']}
            del ds[variableName].attrs['_FillValue']
        else:
            encodingDict[variableName] = {'_FillValue': None}

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

    ds = xarray.open_dataset(args.inFileName, mask_and_scale=False)
    if 'nVertLevels' in ds.dims:
        ds = ds.rename({'nVertLevels': 'depth'})

        dsCoord = xarray.open_dataset(coordFileName, mask_and_scale=False)
        dsCoord = dsCoord.rename({'nVertLevels': 'depth'})

        ds.coords['depth'] = ('depth',
                              compute_depth(dsCoord.refBottomDepth))
        ds.depth.attrs['units'] = 'meters'
        ds.depth.attrs['positive'] = 'down'
        ds.depth.attrs['standard_name'] = 'depth'

        ds.depth.attrs['long_name'] = 'reference depth of the center of ' \
                                      'each vertical level'

        for varName in ds.data_vars:
            var = ds[varName]
            if 'depth' in var.dims:
                var = var.assign_coords(depth=ds.depth)
                ds[varName] = var

    time = datetime.now().strftime('%c')

    history = '{}: {}'.format(time, ' '.join(sys.argv))

    if 'history' in ds.attrs:
        ds.attrs['history'] = '{}\n{}'.format(history,
                                              ds.attrs['history'])
    else:
        ds.attrs['history'] = history

    write_netcdf(ds, args.outFileName)


if __name__ == '__main__':
    main()
