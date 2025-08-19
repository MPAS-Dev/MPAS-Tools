import numpy as np
import xarray as xr

from mpas_tools.ocean.depth import compute_zmid


def compute_vertically_integrated_velocity(
    ds_mesh,
    ds,
    logger=None,
    min_depth=None,
    max_depth=None,
    prefix='timeMonthly_avg_',
    include_bolus=False,
    include_submesoscale=False,
    nedges_chunk=10000,
):
    """
    Compute the vertically integrated normal velocity on mesh edges.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        A dataset containing MPAS mesh variables

    ds : xarray.Dataset
        A dataset containing MPAS output variables ``normalVelocity`` and
        ``layerThickness`` among others, possibly with a ``prefix``

    logger : logging.Logger, optional
        A logger for the output if not stdout

    min_depth : float, optional
        The minimum depth (positive up) to compute BSF over

    max_depth : float, optional
        The maximum depth (positive up) to compute BSF over

    prefix : str, optional
        The prefix on the ``normalVelocity`` and ``layerThickness`` variables

    include_bolus : bool, optional
        Whether to include the GM bolus velocity in the computation

    include_submesoscale : bool, optional
        Whether to include the submesoscale velocity in the computation

    nedges_chunk : int, optional
        The number of edges to chunk the dataset by. This is useful for large
        datasets to avoid memory issues. The default is 10000. Set this to
        ``None`` to disable chunking.

    Returns
    -------
    vert_integ_velocity : xarray.DataArray
        The vertically integrated velocity on the mesh edges
    """
    if nedges_chunk is not None:
        ds = ds.chunk({'nEdges': nedges_chunk})

    if logger:
        logger.info('    Identifying inner edges.')
    cells_on_edge = ds_mesh.cellsOnEdge - 1
    inner_edges = np.logical_and(
        cells_on_edge.isel(TWO=0) >= 0, cells_on_edge.isel(TWO=1) >= 0
    )

    # convert from boolean mask to indices
    inner_edges = np.flatnonzero(inner_edges.values)

    cell0 = cells_on_edge.isel(nEdges=inner_edges, TWO=0)
    cell1 = cells_on_edge.isel(nEdges=inner_edges, TWO=1)
    n_vert_levels = ds.sizes['nVertLevels']

    layer_thickness = ds[f'{prefix}layerThickness']
    max_level_cell = ds_mesh.maxLevelCell - 1
    if 'minLevelCell' in ds_mesh:
        min_level_cell = ds_mesh.minLevelCell - 1
    else:
        min_level_cell = xr.zeros_like(max_level_cell)

    vert_index = xr.DataArray.from_dict(
        {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)}
    )

    if min_depth is not None or max_depth is not None:
        z_mid_edge = _compute_zmid_edge(
            ds, ds_mesh, prefix, logger, cell0, cell1
        )
    else:
        z_mid_edge = None

    if logger:
        logger.info('    Adding up constituents of normal velocity.')
    normal_velocity = ds[f'{prefix}normalVelocity']
    if include_bolus:
        normal_velocity += ds[f'{prefix}normalGMBolusVelocity']
    if include_submesoscale:
        normal_velocity += ds[f'{prefix}normalMLEvelocity']
    normal_velocity = normal_velocity.isel(nEdges=inner_edges)

    thickness0 = layer_thickness.isel(nCells=cell0)
    thickness1 = layer_thickness.isel(nCells=cell1)
    layer_thickness_edge = 0.5 * (thickness0 + thickness1)

    if logger:
        logger.info('    Computing valid cell mask.')

    valid_cells = np.logical_and(
        vert_index >= min_level_cell, vert_index <= max_level_cell
    ).transpose('nCells', 'nVertLevels')
    valid_edges = np.logical_and(
        valid_cells.isel(nCells=cell0), valid_cells.isel(nCells=cell1)
    )

    if logger:
        logger.info('    Applying depth masks.')
    masks = [valid_edges]
    if min_depth is not None:
        masks.append(z_mid_edge <= min_depth)
    if max_depth is not None:
        masks.append(z_mid_edge >= max_depth)
    for mask in masks:
        normal_velocity = normal_velocity.where(mask)
        layer_thickness_edge = layer_thickness_edge.where(mask)

    if logger:
        logger.info('    Computing vertically integrated velocity.')
    vert_integ_velocity = np.zeros(ds_mesh.sizes['nEdges'], dtype=float)
    inner_vert_integ_vel = (layer_thickness_edge * normal_velocity).sum(
        dim='nVertLevels'
    )
    vert_integ_velocity[inner_edges] = inner_vert_integ_vel.values

    if logger:
        logger.info('    Finalizing vertically integrated velocity.')
    vert_integ_velocity = xr.DataArray(vert_integ_velocity, dims=('nEdges',))

    return vert_integ_velocity


def _compute_zmid_edge(ds, ds_mesh, prefix, logger, cell0, cell1):
    layer_thickness = ds[f'{prefix}layerThickness']

    if logger:
        logger.info('    Computing z_mid and z_mid_edge.')
    z_mid = compute_zmid(
        ds_mesh.bottomDepth, ds_mesh.maxLevelCell, layer_thickness
    )
    z_mid_edge = 0.5 * (z_mid.isel(nCells=cell0) + z_mid.isel(nCells=cell1))
    return z_mid_edge
