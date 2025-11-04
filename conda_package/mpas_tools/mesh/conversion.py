import os
import shutil
from tempfile import TemporaryDirectory, mkdtemp

import numpy as np
import xarray as xr

import mpas_tools.io
from mpas_tools.io import write_netcdf
from mpas_tools.logging import check_call


def convert(dsIn, graphInfoFileName=None, logger=None, dir=None):
    """
    Convert an input mesh to a valid MPAS mesh using ``MpasMeshConverter.x``.

    This function ensures the mesh is fully compliant with the
    `MPAS mesh specification <https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf>`_.
    It writes the input dataset to a temporary file, runs the converter, and
    loads the output as an xarray.Dataset.

    Parameters
    ----------
    dsIn : xarray.Dataset
        Input dataset describing the mesh to convert.

    graphInfoFileName : str, optional
        Path to save the generated ``graph.info`` file. If not provided,
        the file is not saved.

    logger : logging.Logger, optional
        Logger for capturing output from the converter.

    dir : str, optional
        Directory in which to create a temporary working directory.

    Returns
    -------
    dsOut : xarray.Dataset
        The converted MPAS mesh dataset.

    Notes
    -----
    - Requires ``MpasMeshConverter.x`` to be available in the system path.
    - Temporary files are created and deleted automatically.
    """  # noqa: E501
    if dir is not None:
        dir = os.path.abspath(dir)

    tempdir = mkdtemp(dir=dir)
    inFileName = f'{tempdir}/mesh_in.nc'
    write_netcdf(dsIn, inFileName)

    outFileName = f'{tempdir}/mesh_out.nc'

    if graphInfoFileName is not None:
        graphInfoFileName = os.path.abspath(graphInfoFileName)

    outDir = os.path.dirname(outFileName)

    check_call(['MpasMeshConverter.x', inFileName, outFileName], logger)

    dsOut = xr.open_dataset(outFileName)
    dsOut.load()

    if graphInfoFileName is not None:
        shutil.copyfile(f'{outDir}/graph.info', graphInfoFileName)

    return dsOut


def cull(
    dsIn,
    dsMask=None,
    dsInverse=None,
    dsPreserve=None,
    graphInfoFileName=None,
    logger=None,
    dir=None,
):
    """
    Cull cells from a mesh using ``MpasCellCuller.x``.

    The function removes cells based on the ``cullCell`` field in the input
    dataset and/or provided mask datasets. Masks are merged as follows:
    - ``cullCell``, ``dsMask``, and ``dsInverse`` are combined (union).
    - ``dsPreserve`` indicates cells that must not be culled.

    Parameters
    ----------
    dsIn : xarray.Dataset
        Input mesh dataset, optionally with a ``cullCell`` field.

    dsMask : xarray.Dataset or list, optional
        Dataset(s) with region masks (1 where cells should be culled).

    dsInverse : xarray.Dataset or list, optional
        Dataset(s) with region masks (0 where cells should be culled).

    dsPreserve : xarray.Dataset or list, optional
        Dataset(s) with region masks (1 where cells should NOT be culled).

    graphInfoFileName : str, optional
        Path to save the generated ``culled_graph.info`` file.

    logger : logging.Logger, optional
        Logger for capturing output from the culler.

    dir : str, optional
        Directory in which to create a temporary working directory.

    Returns
    -------
    dsOut : xarray.Dataset
        The culled mesh dataset.

    Notes
    -----
    - Requires ``MpasCellCuller.x`` to be available in the system path.
    - Temporary files are created and deleted automatically.
    """
    if dir is not None:
        dir = os.path.abspath(dir)
    tempdir = mkdtemp(dir=dir)
    inFileName = f'{tempdir}/ds_in.nc'
    dsIn = _masks_to_int(dsIn)
    write_netcdf(dsIn, inFileName)
    outFileName = f'{tempdir}/ds_out.nc'

    args = ['MpasCellCuller.x', inFileName, outFileName]

    if dsMask is not None:
        if not isinstance(dsMask, list):
            dsMask = [dsMask]
        for index, ds in enumerate(dsMask):
            ds = _masks_to_int(ds)
            fileName = f'{tempdir}/mask{index}.nc'
            write_netcdf(ds, fileName)
            args.extend(['-m', fileName])

    if dsInverse is not None:
        if not isinstance(dsInverse, list):
            dsInverse = [dsInverse]
        for index, ds in enumerate(dsInverse):
            ds = _masks_to_int(ds)
            fileName = f'{tempdir}/inverse{index}.nc'
            write_netcdf(ds, fileName)
            args.extend(['-i', fileName])

    if dsPreserve is not None:
        if not isinstance(dsPreserve, list):
            dsPreserve = [dsPreserve]
        for index, ds in enumerate(dsPreserve):
            ds = _masks_to_int(ds)
            fileName = f'{tempdir}/preserve{index}.nc'
            write_netcdf(ds, fileName)
            args.extend(['-p', fileName])

    if graphInfoFileName is not None:
        graphInfoFileName = os.path.abspath(graphInfoFileName)

    outDir = os.path.dirname(outFileName)

    check_call(args=args, logger=logger)

    dsOut = xr.open_dataset(outFileName)
    dsOut.load()

    if graphInfoFileName is not None:
        shutil.copyfile(f'{outDir}/culled_graph.info', graphInfoFileName)

    return dsOut


def mask(dsMesh, fcMask=None, logger=None, dir=None, cores=1):
    """
    Create region masks on an MPAS mesh using ``compute_mpas_region_masks``.

    Parameters
    ----------
    dsMesh : xarray.Dataset
        MPAS mesh dataset.

    fcMask : geometric_features.FeatureCollection, optional
        Feature collection with regions to mask.

    logger : logging.Logger, optional
        Logger for capturing output.

    dir : str, optional
        Directory in which to create a temporary working directory.

    cores : int, optional
        Number of processes for Python multiprocessing.

    Returns
    -------
    dsMask : xarray.Dataset
        Dataset containing the computed masks.

    Notes
    -----
    - Requires ``compute_mpas_region_masks`` command-line tool.
    - Temporary files are created and deleted automatically.
    """
    if dir is not None:
        dir = os.path.abspath(dir)
    with TemporaryDirectory(dir=dir) as tempdir:
        inFileName = f'{tempdir}/mesh_in.nc'
        write_netcdf(dsMesh, inFileName)
        outFileName = f'{tempdir}/mask_out.nc'

        geojsonFileName = f'{tempdir}/mask.geojson'
        fcMask.to_geojson(geojsonFileName)
        args = [
            'compute_mpas_region_masks',
            '-m',
            inFileName,
            '-o',
            outFileName,
            '-g',
            geojsonFileName,
            '-t',
            'cell',
            '--process_count',
            f'{cores}',
            '--format',
            mpas_tools.io.default_format,
        ]
        if mpas_tools.io.default_engine is not None:
            args.extend(['--engine', mpas_tools.io.default_engine])

        check_call(args=args, logger=logger)

        dsOut = xr.open_dataset(outFileName)
        dsOut.load()

    return dsOut


def _masks_to_int(dsIn):
    """Convert mask variables to int32 type as required by the cell culler."""
    var_list = [
        'regionCellMasks',
        'transectCellMasks',
        'cullCell',
        'cellSeedMask',
    ]
    dsOut = xr.Dataset(dsIn, attrs=dsIn.attrs)
    for var in var_list:
        if var in dsIn:
            dsOut[var] = dsIn[var].astype(np.int32)

    return dsOut
