#!/usr/bin/env python

"""
Add a 3D coordinate "zMid" to an MPAS-Ocean output file that defines the
positive-up vertical location of each cell center.
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


def compute_zmid(bottomDepth, maxLevelCell, layerThickness):
    """
    Computes zMid given data arrays for bottomDepth, maxLevelCell and
    layerThickness

    Parameters
    ----------
    bottomDepth : ``xarray.DataArray``
        the depth of the ocean bottom (positive)

    maxLevelCell : ``xarray.DataArray``
        the 1-based vertical index of the bottom of the ocean

    layerThickness : ``xarray.DataArray``
        the thickness of MPAS-Ocean layers (possibly as a function of time)

    Returns
    -------
    zMid : ``xarray.DataArray``
        the vertical coordinate defining the middle of each layer, masked below
        the bathymetry
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    nDepth = layerThickness.sizes['depth']

    vertIndex = \
        xarray.DataArray.from_dict({'dims': ('depth',),
                                    'data': numpy.arange(nDepth)})

    layerThickness = layerThickness.where(vertIndex < maxLevelCell)

    thicknessSum = layerThickness.sum(dim='depth')
    thicknessCumSum = layerThickness.cumsum(dim='depth')
    zSurface = -bottomDepth+thicknessSum

    zLayerBot = zSurface - thicknessCumSum

    zMid = zLayerBot + 0.5*layerThickness

    zMid = zMid.where(vertIndex < maxLevelCell)
    zMid = zMid.transpose('Time', 'nCells', 'depth')

    return zMid


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--coordFileName", dest="coordFileName",
                        type=str, required=False,
                        help="A MPAS-Ocean file with bottomDepth, maxLevelCell"
                             "and layerThickness but not zMid")
    parser.add_argument("-i", "--inFileName", dest="inFileName", type=str,
                        required=True,
                        help="An input MPAS-Ocean file that zMid should be"
                             "added to, used for coords if another file is"
                             "not provided via -c.")
    parser.add_argument("-o", "--outFileName", dest="outFileName", type=str,
                        required=True,
                        help="An output MPAS-Ocean file with zMid added")
    args = parser.parse_args()

    if args.coordFileName:
        coordFileName = args.coordFileName
    else:
        coordFileName = args.inputFileName

    dsCoord = xarray.open_dataset(coordFileName)
    dsCoord = dsCoord.rename({'nVertLevels': 'depth'})

    ds = xarray.open_dataset(args.inFileName)
    ds = ds.rename({'nVertLevels': 'depth'})

    ds.coords['zMid'] = compute_zmid(dsCoord.bottomDepth, dsCoord.maxLevelCell,
                                     dsCoord.layerThickness)
    ds.zMid.attrs['units'] = 'meters'
    ds.zMid.attrs['positive'] = 'up'
    ds.zMid.attrs['standard_name'] = 'depth'

    for varName in ds.data_vars:
        var = ds[varName]
        if 'nCells' in var.dims and 'depth' in var.dims:
            var = var.assign_coords(zMid=ds.zMid)
            ds[varName] = var

    if 'history' in ds.attrs:
        ds.attrs['history'] = '{}\n{}'.format(' '.join(sys.argv),
                                              ds.attrs['history'])
    else:
        ds.attrs['history'] = ' '.join(sys.argv)

    write_netcdf(ds, args.outFileName)


if __name__ == '__main__':
    main()
