import os
import xarray
from tempfile import TemporaryDirectory, mkdtemp
import shutil

import mpas_tools.io
from mpas_tools.io import write_netcdf
from mpas_tools.logging import check_call


def convert(dsIn, graphInfoFileName=None, logger=None, dir=None):
    """
    Use ``MpasMeshConverter.x`` to convert an input mesh to a valid MPAS
    mesh that is fully compliant with the MPAS mesh specification.
    https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf

    Parameters
    ----------
    dsIn : xarray.Dataset
        A data set to convert

    graphInfoFileName : str, optional
        A file path (relative or absolute) where the graph file (typically
        ``graph.info`` should be written out.  By default, ``graph.info`` is
        not saved.

    logger : logging.Logger, optional
        A logger for the output if not stdout

    dir : str, optional
        A directory in which a temporary directory will be added with files
        produced during conversion and then deleted upon completion.

    Returns
    -------
    dsOut : xarray.Dataset
        The MPAS mesh
    """
    if dir is not None:
        dir = os.path.abspath(dir)

    tempdir = mkdtemp(dir=dir)
    inFileName = '{}/mesh_in.nc'.format(tempdir)
    write_netcdf(dsIn, inFileName)

    outFileName = '{}/mesh_out.nc'.format(tempdir)

    if graphInfoFileName is not None:
        graphInfoFileName = os.path.abspath(graphInfoFileName)

    outDir = os.path.dirname(outFileName)

    check_call(['MpasMeshConverter.x', inFileName, outFileName], logger)

    dsOut = xarray.open_dataset(outFileName)
    dsOut.load()

    if graphInfoFileName is not None:
        shutil.copyfile('{}/graph.info'.format(outDir),
                        graphInfoFileName)

    return dsOut


def cull(dsIn, dsMask=None, dsInverse=None, dsPreserve=None,
         graphInfoFileName=None, logger=None, dir=None):
    """
    Use ``MpasCellCuller.x`` to cull cells from a mesh based on the
    ``cullCell`` field in the input file or DataSet and/or the provided masks.
    ``cullCell``, dsMask and dsInverse are merged together so that the final
    mask is the union of these 3.  The preserve mask is then used to determine
    where cells should *not* be culled.

    Parameters
    ----------
    dsIn : xarray.Dataset
        A data set to cull, possibly with a ``cullCell`` field set to one where
        cells should be removed

    dsMask : xarray.Dataset or list, optional
        A data set (or data sets) with region masks that are 1 where cells
        should be culled

    dsInverse : xarray.Dataset or list, optional
        A data set (or data sets) with region masks that are 0 where cells
        should be culled

    dsPreserve : xarray.Dataset or list, optional
        A data set (or data sets) with region masks that are 1 where cells
        should *not* be culled

    graphInfoFileName : str, optional
        A file path (relative or absolute) where the graph file (typically
        ``culled_graph.info`` should be written out.  By default,
        ``culled_graph.info`` is not saved.

    logger : logging.Logger, optional
        A logger for the output if not stdout

    dir : str, optional
        A directory in which a temporary directory will be added with files
        produced during cell culling and then deleted upon completion.

    Returns
    -------
    dsOut : xarray.Dataset
        The culled mesh

    """
    if dir is not None:
        dir = os.path.abspath(dir)
    tempdir = mkdtemp(dir=dir)
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

    if graphInfoFileName is not None:
        graphInfoFileName = os.path.abspath(graphInfoFileName)

    outDir = os.path.dirname(outFileName)

    check_call(args=args, logger=logger)

    dsOut = xarray.open_dataset(outFileName)
    dsOut.load()

    if graphInfoFileName is not None:
        shutil.copyfile('{}/culled_graph.info'.format(outDir),
                        graphInfoFileName)

    return dsOut


def mask(dsMesh, fcMask=None, logger=None, dir=None, cores=1):
    """
    Use ``compute_mpas_region_masks`` to create a set of region masks either
    from mask feature collections

    Parameters
    ----------
    dsMesh : xarray.Dataset, optional
        An MPAS mesh on which the masks should be created

    fcMask : geometric_features.FeatureCollection, optional
        A feature collection containing features to use to create the mask

    logger : logging.Logger, optional
        A logger for the output if not stdout

    dir : str, optional
        A directory in which a temporary directory will be added with files
        produced during mask creation and then deleted upon completion.

    cores : int, optional
        The number of cores to use for python multiprocessing

    Returns
    -------
    dsMask : xarray.Dataset
        The masks

    """
    if dir is not None:
        dir = os.path.abspath(dir)
    with TemporaryDirectory(dir=dir) as tempdir:
        inFileName = f'{tempdir}/mesh_in.nc'
        write_netcdf(dsMesh, inFileName)
        outFileName = f'{tempdir}/mask_out.nc'

        geojsonFileName = f'{tempdir}/mask.geojson'
        fcMask.to_geojson(geojsonFileName)
        args = ['compute_mpas_region_masks',
                '-m', inFileName,
                '-o', outFileName,
                '-g', geojsonFileName,
                '-t', 'cell',
                '--process_count', f'{cores}',
                '--format', mpas_tools.io.default_format,
                ]
        if mpas_tools.io.default_engine is not None:
            args.extend(['--engine', mpas_tools.io.default_engine])

        check_call(args=args, logger=logger)

        dsOut = xarray.open_dataset(outFileName)
        dsOut.load()

    return dsOut
