import xarray
import numpy
import argparse
import sys
from datetime import datetime
import netCDF4

from mpas_tools.io import write_netcdf


def compute_depth(refBottomDepth):
    """
    Computes depth and depth bounds given refBottomDepth

    Parameters
    ----------
    refBottomDepth : xarray.DataArray
        the depth of the bottom of each vertical layer in the initial state
        (perfect z-level coordinate)

    Returns
    -------
    depth : numpy.ndarray
        the vertical coordinate defining the middle of each layer

    depth_bnds : numpy.ndarray
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


def compute_zmid(bottomDepth, maxLevelCell, layerThickness,
                 depth_dim='nVertLevels'):
    """
    Computes zMid given data arrays for bottomDepth, maxLevelCell and
    layerThickness

    Parameters
    ----------
    bottomDepth : xarray.DataArray
        the depth of the ocean bottom (positive)

    maxLevelCell : xarray.DataArray
        the 1-based vertical index of the bottom of the ocean

    layerThickness : xarray.DataArray
        the thickness of MPAS-Ocean layers (possibly as a function of time)

    depth_dim : str, optional
        the name of the vertical dimension

    Returns
    -------
    zMid : xarray.DataArray
        the vertical coordinate defining the middle of each layer, masked below
        the bathymetry
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    nDepth = layerThickness.sizes[depth_dim]

    vertIndex = \
        xarray.DataArray.from_dict({'dims': (depth_dim,),
                                    'data': numpy.arange(nDepth)})

    layerThickness = layerThickness.where(vertIndex < maxLevelCell)

    thicknessSum = layerThickness.sum(dim=depth_dim)
    thicknessCumSum = layerThickness.cumsum(dim=depth_dim)
    zSurface = -bottomDepth+thicknessSum

    zLayerBot = zSurface - thicknessCumSum

    zMid = zLayerBot + 0.5*layerThickness

    zMid = zMid.where(vertIndex < maxLevelCell)
    zMid = zMid.transpose('Time', 'nCells', depth_dim)

    return zMid


def add_depth(inFileName, outFileName, coordFileName=None):
    """
    Add a 1D depth coordinate to an MPAS-Ocean file.

    Parameters
    ----------
    inFileName : str
        An input MPAS-Ocean file that depth should be added to, used for coords
        if another file is not provided via ``coordFileName``.

    outFileName : str
        An output MPAS-Ocean file with depth added

    coordFileName : str, optional
        An MPAS-Ocean file with ``refBottomDepth``
    """

    if coordFileName is not None:
        coordFileName = inFileName

    ds = xarray.open_dataset(inFileName, mask_and_scale=False)
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

    write_netcdf(ds, outFileName)


def main_add_depth():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--coordFileName", dest="coordFileName",
                        type=str, required=False,
                        help="An MPAS-Ocean file with refBottomDepth")
    parser.add_argument("-i", "--inFileName", dest="inFileName", type=str,
                        required=True,
                        help="An input MPAS-Ocean file that depth should be"
                             "added to, used for coords if another file is"
                             "not provided via -c.")
    parser.add_argument("-o", "--outFileName", dest="outFileName", type=str,
                        required=True,
                        help="An output MPAS-Ocean file with depth added")
    args = parser.parse_args()

    add_depth(args.inFileName, args.outFileName,
              coordFileName=args.coordFileName)


def add_zmid(inFileName, outFileName, coordFileName=None):
    """
    Add a 3D, time-independent depth coordinate to an MPAS-Ocean file.

    Parameters
    ----------
    inFileName : str
        An input MPAS-Ocean file that ``zMid`` should be added to, used for
        coords if another file is not provided via ``coordFileName``.

    outFileName : str
        An output MPAS-Ocean file with ``zMid`` added

    coordFileName : str, optional
        An MPAS-Ocean file with ``bottomDepth``, ``maxLevelCell`` and
        ``layerThickness`` but not ``zMid``
    """
    if coordFileName is None:
        coordFileName = inFileName

    ds = xarray.open_dataset(inFileName, mask_and_scale=False)
    if 'nVertLevels' in ds.dims:
        ds = ds.rename({'nVertLevels': 'depth'})

        # dsCoord doesn't have masking disabled because we want it for zMid
        dsCoord = xarray.open_dataset(coordFileName)
        dsCoord = dsCoord.rename({'nVertLevels': 'depth'})

        ds.coords['zMid'] = compute_zmid(dsCoord.bottomDepth,
                                         dsCoord.maxLevelCell,
                                         dsCoord.layerThickness,
                                         depth_dim='depth')
        fillValue = netCDF4.default_fillvals['f8']
        ds.coords['zMid'] = ds.zMid.where(ds.zMid.notnull(), other=fillValue)
        ds.zMid.attrs['units'] = 'meters'
        ds.zMid.attrs['positive'] = 'up'
        ds.zMid.attrs['_FillValue'] = fillValue

        for varName in ds.data_vars:
            var = ds[varName]
            if 'nCells' in var.dims and 'depth' in var.dims:
                var = var.assign_coords(zMid=ds.zMid)
                ds[varName] = var

    time = datetime.now().strftime('%c')

    history = '{}: {}'.format(time, ' '.join(sys.argv))

    if 'history' in ds.attrs:
        ds.attrs['history'] = '{}\n{}'.format(history,
                                              ds.attrs['history'])
    else:
        ds.attrs['history'] = history

    write_netcdf(ds, outFileName)


def main_add_zmid():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--coordFileName", dest="coordFileName",
                        type=str, required=False,
                        help="An MPAS-Ocean file with bottomDepth, maxLevelCell"
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

    add_zmid(args.inFileName, args.outFileName,
             coordFileName=args.coordFileName)


def write_time_varying_zmid(inFileName, outFileName, coordFileName=None,
                            prefix=None):
    """
    Add a 3D, time-independent depth coordinate to an MPAS-Ocean file.

    Parameters
    ----------
    inFileName : str
        An input MPAS-Ocean file with some form of ``layerThickness``, and also
        ``bottomDepth`` and ``maxLevelCell`` if no ``coordFileName``
        is provided.

    outFileName : str
        An output MPAS-Ocean file with ``zMid`` for each time in the input file

    coordFileName : str, optional
        An MPAS-Ocean file with ``bottomDepth`` and ``maxLevelCell``

    prefix : str, optional
        A prefix on ``layerThickness`` (in) and ``zMid`` (out), such as
        ``timeMonthly_avg_``
    """
    if coordFileName is None:
        coordFileName = inFileName

    dsCoord = xarray.open_dataset(coordFileName)
    dsCoord = dsCoord.rename({'nVertLevels': 'depth'})

    dsIn = xarray.open_dataset(inFileName)
    dsIn = dsIn.rename({'nVertLevels': 'depth'})
    inVarName = '{}layerThickness'.format(prefix)
    outVarName = '{}zMid'.format(prefix)
    layerThickness = dsIn[inVarName]

    zMid = compute_zmid(dsCoord.bottomDepth, dsCoord.maxLevelCell,
                        layerThickness, depth_dim='depth')

    dsOut = xarray.Dataset()
    dsOut[outVarName] = zMid
    fillValue = netCDF4.default_fillvals['f8']
    dsOut[outVarName] = dsOut[outVarName].where(dsOut[outVarName].notnull(),
                                                other=fillValue)
    dsOut[outVarName].attrs['units'] = 'meters'
    dsOut[outVarName].attrs['positive'] = 'up'
    dsOut[outVarName].attrs['_FillValue'] = fillValue

    time = datetime.now().strftime('%c')

    history = '{}: {}'.format(time, ' '.join(sys.argv))
    dsOut.attrs['history'] = history

    write_netcdf(dsOut, outFileName)


def main_write_time_varying_zmid():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--coordFileName", dest="coordFileName",
                        type=str, required=False,
                        help="An MPAS-Ocean file with bottomDepth and "
                             "maxLevelCell")
    parser.add_argument("-i", "--inFileName", dest="inFileName", type=str,
                        required=True,
                        help="An input MPAS-Ocean file with some form of"
                             "layerThickness, and also bottomDepth and"
                             "maxLevelCell if no coordinate file is provided.")
    parser.add_argument("-o", "--outFileName", dest="outFileName", type=str,
                        required=True,
                        help="An output MPAS-Ocean file with zMid for each"
                             "time in the input file")
    parser.add_argument("-p", "--prefix", dest="prefix", type=str,
                        required=False, default="",
                        help="A prefix on layerThickness (in) and zMid (out),"
                             "such as 'timeMonthly_avg_'")
    args = parser.parse_args()

    write_time_varying_zmid(
        args.inFileName, args.outFileName, coordFileName=args.coordFileName,
        prefix=args.prefix)
