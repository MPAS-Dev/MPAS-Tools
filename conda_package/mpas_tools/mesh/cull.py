import numpy as np
import xarray as xr
from scipy.spatial import KDTree

from mpas_tools.io import write_netcdf


def write_map_culled_to_base(base_mesh_filename, culled_mesh_filename,
                             out_filename, workers=-1):
    """
    Write out a file with maps from cells, edges and vertices on a culled mesh
    to the same elements on a base mesh.  All elements in the culled mesh must
    be in the base mesh.

    Parameters
    ----------
    base_mesh_filename : str
        A file with the horizontal MPAS mesh before culling

    culled_mesh_filename : str
        A file with the culled horizonal MPAS mesh

    out_filename : str
        A file to which the maps should be written.  The dataset will include
        ``mapCulledToBaseCell``, ``mapCulledToBaseEdge`` and
        ``mapCulledToBaseVertex``, each of which contains the base mesh index
        corresponding to the same element from the culled mesh.

    workers : int, optional
        The number of threads to use to query base mesh elements. The default
        is all available threads (``workers=-1``)
    """
    ds_base = xr.open_dataset(base_mesh_filename)
    ds_culled = xr.open_dataset(culled_mesh_filename)
    ds_map_culled_to_base = map_culled_to_base(ds_base, ds_culled,
                                               workers=workers)
    write_netcdf(ds_map_culled_to_base, out_filename)


def map_culled_to_base(ds_base, ds_culled, workers=-1):
    """
    Create maps from cells, edges and vertices on a culled mesh to the same
    elements on a base mesh.  All elements in the culled mesh must be in the
    base mesh.

    Parameters
    ----------
    ds_base : xarray.Dataset
        The horizontal MPAS mesh before culling

    ds_culled : xarray.Dataset
        The culled horizonal MPAS mesh

    workers : int, optional
        The number of threads to use to query base mesh elements. The default
        is all available threads (``workers=-1``)

    Returns
    -------
    ds_map_culled_to_base : xarray.Dataset
        A dataset with ``mapCulledToBaseCell``, ``mapCulledToBaseEdge`` and
        ``mapCulledToBaseVertex``, each of which contains the base mesh index
        corresponding to the same element from the culled mesh.
    """
    ds_map_culled_to_base = xr.Dataset()
    for dim, suffix in [('nCells', 'Cell'),
                        ('nEdges', 'Edge'),
                        ('nVertices', 'Vertex')]:
        _map_culled_to_base_grid_type(ds_base, ds_culled,
                                      ds_map_culled_to_base, dim, suffix,
                                      workers)

    return ds_map_culled_to_base


def write_culled_dataset(in_filename, out_filename, base_mesh_filename,
                         culled_mesh_filename,
                         map_culled_to_base_filename=None, workers=-1,
                         logger=None):
    """
    Create a new dataset in which all fields from ``ds`` have been culled
    from the base mesh to the culled mesh.  Fields present in
    ``ds_culled_mesh`` are copied over rather than culled from ``ds``.

    Parameters
    ----------
    in_filename : str
        A file containing an MPAS dataset to cull

    output_filename : str
        A file to write the culled MPAS dataset to

    base_mesh_filename : str
        A file with the horizontal MPAS mesh before culling

    culled_mesh_filename : str
        A file with the culled horizonal MPAS mesh

    map_culled_to_base_filename : str, optional
        A file with an existing map from the base to the culled mesh created
        with ``write_map_culled_to_base()`` or ``map_culled_to_base()``. The
        dataset will be created (but not returned or saved to disk) if it is
        not passed as an argument.

    workers : int, optional
        The number of threads to use to query base mesh elements. The default
        is all available threads (``workers=-1``)

    logger : logging.Logger, optional
        A logger for the output
    """
    ds = xr.open_dataset(in_filename)
    ds_base_mesh = xr.open_dataset(base_mesh_filename)
    ds_culled_mesh = xr.open_dataset(culled_mesh_filename)
    if map_culled_to_base_filename is None:
        ds_map_culled_to_base = None
    else:
        ds_map_culled_to_base = xr.open_dataset(map_culled_to_base_filename)

    ds_culled = cull_dataset(
        ds=ds, ds_base_mesh=ds_base_mesh, ds_culled_mesh=ds_culled_mesh,
        ds_map_culled_to_base=ds_map_culled_to_base,
        workers=workers, logger=logger)
    write_netcdf(ds_culled, out_filename)


def cull_dataset(ds, ds_base_mesh, ds_culled_mesh, ds_map_culled_to_base=None,
                 workers=-1, logger=None):
    """
    Create a new dataset in which all fields from ``ds`` have been culled
    from the base mesh to the culled mesh.  Fields present in
    ``ds_culled_mesh`` are copied over rather than culled from ``ds``.

    Parameters
    ----------
    ds : xarray.Dataset
        An MPAS dataset to cull

    ds_base_mesh : xarray.Dataset
        The horizontal MPAS mesh before culling

    ds_culled_mesh : xarray.Dataset
        The culled horizonal MPAS mesh

    ds_map_culled_to_base : xarray.Dataset, optional
        An existing map from the base to the culled mesh created with
        ``write_map_culled_to_base()`` or ``map_culled_to_base()``. The dataset
        will be created (but not returned or saved to disk) if it is not passed
        as an argument.

    workers : int, optional
        The number of threads to use to query base mesh elements. The default
        is all available threads (``workers=-1``)

    logger : logging.Logger, optional
        A logger for the output

    Returns
    -------
    ds_culled : xarray.Dataset
        An culled MPAS dataset
    """
    if ds_map_culled_to_base is None:
        if logger is not None:
            logger.info('Creating culled-to-base mapping')
        ds_map_culled_to_base = map_culled_to_base(
            ds_base=ds_base_mesh, ds_culled=ds_culled_mesh, workers=workers)

    if logger is not None:
        logger.info('Culling dataset')
    ds_culled = ds
    if 'nCells' in ds_culled.dims:
        ds_culled = ds_culled.isel(
            nCells=ds_map_culled_to_base['mapCulledToBaseCell'].values)
    if 'nEdges' in ds_culled.dims:
        ds_culled = ds_culled.isel(
            nEdges=ds_map_culled_to_base['mapCulledToBaseEdge'].values)
    if 'nVertices' in ds_culled.dims:
        ds_culled = ds_culled.isel(
            nVertices=ds_map_culled_to_base['mapCulledToBaseVertex'].values)

    if logger is not None:
        logger.info('Replacing variables from culled mesh')
    for var in ds.data_vars:
        if var in ds_culled_mesh:
            if logger is not None:
                logger.info(f'  replacing: {var}')
            # replace this field with the version from the culled mesh
            ds_culled[var] = ds_culled_mesh[var]
        else:
            if logger is not None:
                logger.info(f'  keeping:   {var}')

    return ds_culled


def _map_culled_to_base_grid_type(ds_base, ds_culled, ds_map_culled_to_base,
                                  dim, suffix, workers):
    x_base = ds_base[f'x{suffix}'].values
    y_base = ds_base[f'y{suffix}'].values
    z_base = ds_base[f'z{suffix}'].values

    x_culled = ds_culled[f'x{suffix}'].values
    y_culled = ds_culled[f'y{suffix}'].values
    z_culled = ds_culled[f'z{suffix}'].values

    # create a map from lat-lon pairs to base-mesh cell indices
    points = np.vstack((x_base, y_base, z_base)).T

    tree = KDTree(points)

    points = np.vstack((x_culled, y_culled, z_culled)).T

    _, culled_to_base_map = tree.query(points, workers=workers)

    ds_map_culled_to_base[f'mapCulledToBase{suffix}'] = \
        ((dim,), culled_to_base_map)
