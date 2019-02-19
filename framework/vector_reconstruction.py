#!/usr/bin/env python

"""
Extract Cartesian (X, Y, Z), zonal and meridional components of an MPAS vector
field, given the field on edge normals.

This tool requires that the field 'coeffs_reconstruct' has been saved to a
NetCDF file.  The simplest way to do this is to include the following stream
in a forward run:

<stream name="vector_reconstruction"
        clobber_mode="truncate"
        type="output"
        output_interval="0000-00-00_00:00:01"
        filename_template="vector_reconstruction.nc">

    <var name="coeffs_reconstruct"/>
</stream>

and run the model for one time step.

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
from datetime import datetime
from dask.diagnostics import ProgressBar


def write_netcdf(ds, fileName, fillValues=netCDF4.default_fillvals):
    encodingDict = {}
    variableNames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variableName in variableNames:
        dtype = ds[variableName].dtype
        for fillType in fillValues:
            if dtype == numpy.dtype(fillType):
                encodingDict[variableName] = \
                    {'_FillValue': fillValues[fillType]}
                break

    delayed_obj = ds.to_netcdf(fileName, encoding=encodingDict, compute=False)

    print('Writing {}'.format(fileName))
    with ProgressBar():
        delayed_obj.compute()


def reconstruct_variable(outVarName, variableOnEdges, dsMesh,
                         coeffs_reconstruct, dsOut, chunkSize=32768):
    nCells = dsMesh.sizes['nCells']
    # nEdgesOnCell = dsMesh.nEdgesOnCell.values
    edgesOnCell = dsMesh.edgesOnCell - 1

    variableOnEdges.load()
    edgesOnCell.load()
    coeffs_reconstruct.load()

    dims = []
    sizes = []
    varIndices = {}
    for dim in variableOnEdges.dims:
        size = variableOnEdges.sizes[dim]
        varIndices[dim] = numpy.arange(size)
        if dim == 'nEdges':
            dim = 'nCells'
            size = nCells
            varIndices['nEdges'] = edgesOnCell
        dims.append(dim)
        sizes.append(size)

    coeffs_reconstruct = coeffs_reconstruct.chunk({'nCells': chunkSize})

    variable = variableOnEdges[varIndices].chunk({'nCells': chunkSize})
    print('Computing {} at edgesOnCell:'.format(outVarName))
    with ProgressBar():
        variable.compute()

    varCart = []

    print('Computing Cartesian conponents:')
    for index, component in enumerate(['X', 'Y', 'Z']):
        var = (coeffs_reconstruct.isel(R3=index)*variable).sum(
            dim='maxEdges').transpose(*dims)
        outName = '{}{}'.format(outVarName, component)
        print(outName)
        with ProgressBar():
            var.compute()
        dsOut[outName] = var
        varCart.append(var)

    latCell = dsMesh.latCell
    lonCell = dsMesh.lonCell
    latCell.load()
    lonCell.load()

    clat = numpy.cos(latCell)
    slat = numpy.sin(latCell)
    clon = numpy.cos(lonCell)
    slon = numpy.sin(lonCell)

    print('Computing zonal and meridional components:')

    outName = '{}Zonal'.format(outVarName)
    zonal = -varCart[0]*slon + varCart[1]*clon
    print(outName)
    with ProgressBar():
        zonal.compute()
    dsOut[outName] = zonal

    outName = '{}Meridional'.format(outVarName)
    merid = -(varCart[0]*clon + varCart[1]*slon)*slat + varCart[2]*clat
    print(outName)
    with ProgressBar():
        merid.compute()
    dsOut[outName] = merid


def main():

    # client = Client(n_workers=1, threads_per_worker=4, memory_limit='10GB')
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-m", "--meshFileName", dest="meshFileName",
                        type=str, required=False,
                        help="An MPAS file with mesh data (edgesOnCell, etc.)")
    parser.add_argument("-w", "--weightsFileName", dest="weightsFileName",
                        type=str, required=False,
                        help="An MPAS file with coeffs_reconstruct ")
    parser.add_argument("-i", "--inFileName", dest="inFileName", type=str,
                        required=True,
                        help="An MPAS file with one or more fields on edges "
                             "to be reconstructed at cell centers.  Used for "
                             "mesh data and/or weights if a separate files "
                             "are not provided.")
    parser.add_argument("-v", "--variables", dest="variables", type=str,
                        required=True,
                        help="A comma-separated list of variables on edges to "
                             "reconstruct")
    parser.add_argument("--outVariables", dest="outVariables", type=str,
                        required=False,
                        help="A comma-separated list of prefixes for output "
                             "variable names")
    parser.add_argument("-o", "--outFileName", dest="outFileName", type=str,
                        required=True,
                        help="An output MPAS file with the reconstructed "
                             "X, Y, Z, zonal and meridional fields")
    args = parser.parse_args()

    if args.meshFileName:
        meshFileName = args.meshFileName
    else:
        meshFileName = args.inFileName

    if args.weightsFileName:
        weightsFileName = args.weightsFileName
    else:
        weightsFileName = args.inFileName

    variables = args.variables.split(',')
    if args.outVariables:
        outVariables = args.outVariables.split(',')
    else:
        outVariables = variables

    dsIn = xarray.open_dataset(args.inFileName, mask_and_scale=False)
    dsMesh = xarray.open_dataset(meshFileName)
    dsWeights = xarray.open_dataset(weightsFileName)
    coeffs_reconstruct = dsWeights.coeffs_reconstruct
    dsOut = xarray.Dataset()

    for inVarName, outVarName in zip(variables, outVariables):
        reconstruct_variable(outVarName, dsIn[inVarName], dsMesh,
                             coeffs_reconstruct, dsOut)

    for attrName in dsIn.attrs:
        dsOut.attrs[attrName] = dsIn.attrs[attrName]

    time = datetime.now().strftime('%c')

    history = '{}: {}'.format(time, ' '.join(sys.argv))

    if 'history' in dsOut.attrs:
        dsOut.attrs['history'] = '{}\n{}'.format(history,
                                                 dsOut.attrs['history'])
    else:
        dsOut.attrs['history'] = history

    write_netcdf(dsOut, args.outFileName)


if __name__ == '__main__':
    main()
