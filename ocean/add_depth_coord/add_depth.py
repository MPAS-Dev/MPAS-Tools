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
    Computes depth and depth bounds given refBottomDepth

    Parameters
    ----------
    refBottomDepth : ``xarray.DataArray``
        the depth of the bottom of each vertical layer in the initial state
        (perfect z-level coordinate)

    Returns
    -------
    depth : ``xarray.DataArray``
        the vertical coordinate defining the middle of each layer
    depth_bnds : ``xarray.DataArray``
        the vertical coordinate defining the top and bottom of each layer
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    refBottomDepth = refBottomDepth.values

    depth_bnds = numpy.zeros((len(refBottomDepth), 2))

    depth_bnds[0, 0] = 0.
    depth_bnds[1:, 0] = refBottomDepth[0:-1]
    depth_bnds[:, 1] = refBottomDepth
    depth = 0.5*(depth_bnds[:, 0] + depth_bnds[:, 1])

    return depth, depth_bnds


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

        depth, depth_bnds = compute_depth(dsCoord.refBottomDepth)
        ds.coords['depth'] = ('depth', depth)
        ds.depth.attrs['long_name'] = 'reference depth of the center of ' \
                                      'each vertical level'
        ds.depth.attrs['standard_name'] = 'depth'
        ds.depth.attrs['units'] = 'meters'
        ds.depth.attrs['axis'] = 'Z'
        ds.depth.attrs['positive'] = 'down'
        ds.depth.attrs['valid_min'] = depth_bnds[0, 0]
        ds.depth.attrs['valid_max'] = depth_bnds[-1, 1]
        ds.depth.attrs['bounds'] = 'depth_bnds'

        ds.coords['depth_bnds'] = (('depth', 'nbnd'), depth_bnds)
        ds.depth_bnds.attrs['long_name'] = 'Gridcell depth interfaces'

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
