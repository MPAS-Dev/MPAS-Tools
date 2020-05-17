from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray
import subprocess
from backports.tempfile import TemporaryDirectory
import shutil

from mpas_tools.io import write_netcdf


def convert(dsIn, graphInfoFileName=None, logger=None):
    '''
    Use ``MpasMeshConverter.x`` to convert an input mesh to a valid MPAS
    mesh that is fully compliant with the MPAS mesh specification.
    https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf

    Parameters
    ----------
    dsIn : ``xarray.Dataset``
        A data set to convert

    graphInfoFileName : str, optional
        A file path (relative or absolute) where the graph file (typically
        ``graph.info`` should be written out.  By default, ``graph.info`` is
        not saved.

    logger : ``logging.Logger``, optional
        A logger for the output if not stdout

    Returns
    -------
    dsOut : ``xarray.Dataset``
        The MPAS mesh
    '''

    with TemporaryDirectory() as tempdir:
        inFileName = '{}/mesh_in.nc'.format(tempdir)
        write_netcdf(dsIn, inFileName)

        outFileName = '{}/mesh_out.nc'.format(tempdir)

        if graphInfoFileName is not None:
            graphInfoFileName = os.path.abspath(graphInfoFileName)

        # go into the directory of the output file so the graph.info file ends
        # up in the same place
        owd = os.getcwd()
        outDir = os.path.dirname(outFileName)
        os.chdir(outDir)
        _call_subprocess(['MpasMeshConverter.x', inFileName, outFileName],
                         logger)
        os.chdir(owd)

        dsOut = xarray.open_dataset(outFileName)
        dsOut.load()

        if graphInfoFileName is not None:
            shutil.copyfile('{}/graph.info'.format(outDir),
                            graphInfoFileName)

    return dsOut


def cull(dsIn, dsMask=None, dsInverse=None, dsPreserve=None,
         graphInfoFileName=None, logger=None):
    '''
    Use ``MpasCellCuller.x`` to cull cells from a mesh based on the
    ``cullCell`` field in the input file or DataSet and/or the provided masks.
    ``cullCell``, dsMask and dsInverse are merged together so that the final
    mask is the union of these 3.  The preserve mask is then used to determine
    where cells should *not* be culled.

    Parameters
    ----------
    dsIn : ``xarray.Dataset``
        A data set to cull, possibly with a ``cullCell`` field set to one where
        cells should be removed

    dsMask : ``xarray.Dataset`` or list, optional
        A data set (or data sets) with region masks that are 1 where cells
        should be culled

    dsInverse : ``xarray.Dataset`` or list, optional
        A data set (or data sets) with region masks that are 0 where cells
        should be culled

    dsPreserve : ``xarray.Dataset`` or list, optional
        A data set (or data sets) with region masks that are 1 where cells
        should *not* be culled

    graphInfoFileName : str, optional
        A file path (relative or absolute) where the graph file (typically
        ``culled_graph.info`` should be written out.  By default,
        ``culled_graph.info`` is not saved.

    logger : ``logging.Logger``, optional
        A logger for the output if not stdout

    Returns
    -------
    dsOut : ``xarray.Dataset``
        The culled mesh

    '''

    with TemporaryDirectory() as tempdir:
        inFileName = '{}/ds_in.nc'.format(tempdir)
        write_netcdf(dsIn, inFileName)
        outFileName = '{}/ds_out.nc'.format(tempdir)

        args = ['MpasCellCuller.x', inFileName, outFileName]

        if dsMask is not None:
            if not isinstance(dsMask, list):
                dsMask = [dsMask]
            for index, ds in enumerate(dsMask):
                fileName = '{}/mask{}.nc'.format(tempdir, index)
                write_netcdf(ds, fileName)
                args.extend(['-m', fileName])

        if dsInverse is not None:
            if not isinstance(dsInverse, list):
                dsInverse = [dsInverse]
            for index, ds in enumerate(dsInverse):
                fileName = '{}/inverse{}.nc'.format(tempdir, index)
                write_netcdf(ds, fileName)
                args.extend(['-i', fileName])

        if dsPreserve is not None:
            if not isinstance(dsPreserve, list):
                dsPreserve = [dsPreserve]
            for index, ds in enumerate(dsPreserve):
                fileName = '{}/preserve{}.nc'.format(tempdir, index)
                write_netcdf(ds, fileName)
                args.extend(['-p', fileName])

        # go into the directory of the output file so the graph.info file ends
        # up in the same place

        if graphInfoFileName is not None:
            graphInfoFileName = os.path.abspath(graphInfoFileName)

        owd = os.getcwd()
        outDir = os.path.dirname(outFileName)
        os.chdir(outDir)
        _call_subprocess(args, logger)
        os.chdir(owd)

        dsOut = xarray.open_dataset(outFileName)
        dsOut.load()

        if graphInfoFileName is not None:
            shutil.copyfile('{}/culled_graph.info'.format(outDir),
                            graphInfoFileName)

    return dsOut


def mask(dsMesh, fcMask=None, fcSeed=None, positiveLon=False, logger=None):
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

    logger : ``logging.Logger``, optional
        A logger for the output if not stdout

    Returns
    -------
    dsMask : ``xarray.Dataset``
        The masks

    '''

    with TemporaryDirectory() as tempdir:
        inFileName = '{}/mesh_in.nc'.format(tempdir)
        write_netcdf(dsMesh, inFileName)
        outFileName = '{}/mesh_out.nc'.format(tempdir)

        args = ['MpasMaskCreator.x', inFileName, outFileName]

        if fcMask is not None:
            fileName = '{}/mask.geojson'.format(tempdir)
            fcMask.to_geojson(fileName)
            args.extend(['-f', fileName])

        if fcSeed is not None:
            fileName = '{}/seed.geojson'.format(tempdir)
            fcSeed.to_geojson(fileName)
            args.extend(['-s', fileName])

        if positiveLon:
            args.append('--positive_lon')

        _call_subprocess(args, logger)

        dsOut = xarray.open_dataset(outFileName)
        dsOut.load()

    return dsOut


def _call_subprocess(args, logger):
    """Call the given subprocess and send the output to the logger"""
    if logger is None:
        subprocess.check_call(args)
    else:
        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if stdout:
            stdout = stdout.decode('utf-8')
            for line in stdout.split('\n'):
                logger.info(line)
        if stderr:
            stderr = stderr.decode('utf-8')
            for line in stderr.split('\n'):
                logger.error(line)

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode,
                                                ' '.join(args))
