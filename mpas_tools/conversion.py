from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray
import subprocess
import tempfile
import shutil

from mpas_tools.io import write_netcdf


def convert(dsIn):
    '''
    Use ``MpasMeshConverter.x`` to convert an input mesh to a valid MPAS
    mesh that is fully compliant with the MPAS mesh specification.
    https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf

    Parameters
    ----------
    dsIn : ``xarray.Dataset``
        A data set to convert

    Returns
    -------
    dsOut : ``xarray.Dataset``
        The MPAS mesh
    '''

    tempFiles = []
    inFileName = _get_temp_path(tempFiles)
    write_netcdf(dsIn, inFileName)

    outFileName = _get_temp_path(tempFiles)

    # go into the directory of the output file so the graph.info file ends
    # up in the same place
    owd = os.getcwd()
    os.chdir(os.path.dirname(outFileName))
    subprocess.check_call(['MpasMeshConverter.x', inFileName, outFileName])
    os.chdir(owd)

    dsOut = xarray.open_dataset(outFileName)
    dsOut.load()
    _remove_temp_files(tempFiles)

    return dsOut


def cull(dsIn, dsMask=None, dsInverse=None, dsPreserve=None,
         graphInfoPath=None):
    '''
    Use ``MpasCellCuller.x`` to cull cells from a mesh based on the
    ``cullCell`` field in the input file or DataSet and/or the provided masks.
    ``cullCell``, dsMask and dsInverse are merged together so that the final
    mask is the union of these 3.  The preserve mask is then used to determine
    where cells should *not* be culled.

    Parameters
    ----------
    dsIn : ``xarray.Dataset``, optional
        A data set to cull, possibly with a ``cullCell`` field set to one where
        cells should be removed

    dsMask : ``xarray.Dataset``, optional
        A data set with region masks that are 1 where cells should be culled

    dsInverse : ``xarray.Dataset``, optional
        A data set with region masks that are 0 where cells should be culled

    dsPreserve : ``xarray.Dataset``, optional
        A data set with region masks that are 1 where cells should *not* be
        culled

    graphInfoPath : str, optional
        A path where the file ``graph.info`` should be written out.  By
        default, ``graph.info`` is written to a temp directory that is deleted.

    Returns
    -------
    dsOut : ``xarray.Dataset``
        The culled mesh

    '''

    tempFiles = []
    inFileName = _get_temp_path(tempFiles)
    write_netcdf(dsIn, inFileName)
    outFileName = _get_temp_path(tempFiles)

    args = ['MpasCellCuller.x', inFileName, outFileName]

    if dsMask is not None:
        fileName = _get_temp_path(tempFiles)
        write_netcdf(dsMask, fileName)
        args.extend(['-m', fileName])

    if dsInverse is not None:
        fileName = _get_temp_path(tempFiles)
        write_netcdf(dsInverse, fileName)
        args.extend(['-i', fileName])

    if dsPreserve is not None:
        fileName = _get_temp_path(tempFiles)
        write_netcdf(dsPreserve, fileName)
        args.extend(['-p', fileName])

    # go into the directory of the output file so the graph.info file ends
    # up in the same place

    if graphInfoPath is not None:
        graphInfoPath = os.path.abspath(graphInfoPath)

    owd = os.getcwd()
    outDir = os.path.dirname(outFileName)
    os.chdir(outDir)
    subprocess.check_call(args)
    os.chdir(owd)

    dsOut = xarray.open_dataset(outFileName)
    dsOut.load()

    if graphInfoPath is not None:
        shutil.copyfile('{}/graph.info'.format(outDir),
                        '{}/graph.info'.format(graphInfoPath))
    _remove_temp_files(tempFiles)

    return dsOut


def mask(dsMesh, fcMask=None, fcSeed=None, positiveLon=False):
    '''
    Use ``MpasMaskCreator.x`` to create a set of region masks either from
    mask feature collecitons or from seed points to be used to flood fill

    Parameters
    ----------
    dsMesh : ``xarray.Dataset``, optional
        An MPAS mesh on which the masks should be created

    fcMask : ``geometric_features.FeatureCollection``, optional
        A feature collection containing features to use to create the mask

    fcSeed : ``geometric_features.FeatureCollection``, optional
        A feature collection with points to use a seeds for a flood fill that
        will create a mask of all cells connected to the seed points

    Returns
    -------
    dsMask : ``xarray.Dataset``
        The masks

    '''

    tempFiles = []
    inFileName = _get_temp_path(tempFiles)
    write_netcdf(dsMesh, inFileName)
    outFileName = _get_temp_path(tempFiles)

    args = ['MpasMaskCreator.x', inFileName, outFileName]

    if fcMask is not None:
        fileName = _get_temp_path(tempFiles, ext='geojson')
        fcMask.to_geojson(fileName)
        args.extend(['-f', fileName])

    if fcSeed is not None:
        fileName = _get_temp_path(tempFiles, ext='geojson')
        fcSeed.to_geojson(fileName)
        args.extend(['-s', fileName])

    if positiveLon:
        args.append('--positive_lon')

    # go into the directory of the output file so the graph.info file ends
    # up in the same place
    owd = os.getcwd()
    os.chdir(os.path.dirname(outFileName))
    subprocess.check_call(args)
    os.chdir(owd)

    dsOut = xarray.open_dataset(outFileName)
    dsOut.load()
    _remove_temp_files(tempFiles)

    return dsOut


def _get_temp_path(tempFiles, ext='nc'):
    '''Returns the name of a temporary NetCDF file'''
    fileName = '{}/{}.{}'.format(tempfile._get_default_tempdir(),
                                 next(tempfile._get_candidate_names()),
                                 ext)
    tempFiles.append(fileName)
    return fileName


def _remove_temp_files(tempFiles):
    for tempFileName in tempFiles:
        os.remove(tempFileName)
