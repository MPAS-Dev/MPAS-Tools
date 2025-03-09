import logging
import sys

import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import xarray as xr

from mpas_tools.ocean.depth import compute_zmid


def compute_barotropic_streamfunction(ds_mesh, ds, logger=None,
                                      min_depth=-5., max_depth=1.e4,
                                      prefix='timeMonthly_avg_',
                                      time_index=0):
    """
    Compute barotropic streamfunction. Returns BSF in Sv on vertices.

    Parameters
    ----------
    ds_mesh : ``xarray.Dataset``
        A dataset containing MPAS mesh variables

    ds : ``xarray.Dataset``
        A dataset containing MPAS output variables ``normalVelocity`` and
        ``layerThickness`` (possibly with a ``prefix``)

    logger : ``logging.Logger``, optional
        A logger for the output if not stdout

    min_depth : float, optional
        The minimum depth (positive down) to compute transport over

    max_depth : float, optional
        The maximum depth (positive down) to compute transport over

    prefix : str, optional
        The prefix on the ``normalVelocity`` and ``layerThickness`` variables

    time_index : int, optional
        The time at which to index ``ds`` (if it has ``Time`` as a dimension)
    """

    useStdout = logger is None
    if useStdout:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.setLevel(logging.INFO)

    inner_edges, transport = _compute_transport(
        ds_mesh, ds, min_depth=min_depth, max_depth=max_depth, prefix=prefix,
        time_index=time_index)
    logger.info('transport computed.')

    nvertices = ds_mesh.sizes['nVertices']

    cells_on_vertex = ds_mesh.cellsOnVertex - 1
    vertices_on_edge = ds_mesh.verticesOnEdge - 1
    is_boundary_cov = cells_on_vertex == -1
    boundary_vertices = np.logical_or(is_boundary_cov.isel(vertexDegree=0),
                                      is_boundary_cov.isel(vertexDegree=1))
    boundary_vertices = np.logical_or(boundary_vertices,
                                      is_boundary_cov.isel(vertexDegree=2))

    # convert from boolean mask to indices
    boundary_vertices = np.flatnonzero(boundary_vertices.values)

    n_boundary_vertices = len(boundary_vertices)
    n_inner_edges = len(inner_edges)

    indices = np.zeros((2, 2 * n_inner_edges + n_boundary_vertices), dtype=int)
    data = np.zeros(2 * n_inner_edges + n_boundary_vertices, dtype=float)

    # The difference between the streamfunction at vertices on an inner
    # edge should be equal to the transport
    v0 = vertices_on_edge.isel(nEdges=inner_edges, TWO=0).values
    v1 = vertices_on_edge.isel(nEdges=inner_edges, TWO=1).values

    ind = np.arange(n_inner_edges)
    indices[0, 2 * ind] = ind
    indices[1, 2 * ind] = v1
    data[2 * ind] = 1.

    indices[0, 2 * ind + 1] = ind
    indices[1, 2 * ind + 1] = v0
    data[2 * ind + 1] = -1.

    # the streamfunction should be zero at all boundary vertices
    ind = np.arange(n_boundary_vertices)
    indices[0, 2 * n_inner_edges + ind] = n_inner_edges + ind
    indices[1, 2 * n_inner_edges + ind] = boundary_vertices
    data[2 * n_inner_edges + ind] = 1.

    rhs = np.zeros(n_inner_edges + n_boundary_vertices, dtype=float)

    # convert to Sv
    ind = np.arange(n_inner_edges)
    rhs[ind] = 1e-6 * transport

    ind = np.arange(n_boundary_vertices)
    rhs[n_inner_edges + ind] = 0.

    matrix = scipy.sparse.csr_matrix(
        (data, indices),
        shape=(n_inner_edges + n_boundary_vertices, nvertices))

    solution = scipy.sparse.linalg.lsqr(matrix, rhs)
    bsf_vertex = xr.DataArray(-solution[0],
                              dims=('nVertices',))

    return bsf_vertex


def _compute_transport(ds_mesh, ds, min_depth, max_depth, prefix,
                       time_index):

    cells_on_edge = ds_mesh.cellsOnEdge - 1
    inner_edges = np.logical_and(cells_on_edge.isel(TWO=0) >= 0,
                                 cells_on_edge.isel(TWO=1) >= 0)

    if 'Time' in ds.dims:
        ds = ds.isel(Time=time_index)

    # convert from boolean mask to indices
    inner_edges = np.flatnonzero(inner_edges.values)

    cell0 = cells_on_edge.isel(nEdges=inner_edges, TWO=0)
    cell1 = cells_on_edge.isel(nEdges=inner_edges, TWO=1)

    normal_velocity = \
        ds[f'{prefix}normalVelocity'].isel(nEdges=inner_edges)
    layer_thickness = ds[f'{prefix}layerThickness']
    layer_thickness_edge = 0.5 * (layer_thickness.isel(nCells=cell0) +
                                  layer_thickness.isel(nCells=cell1))

    n_vert_levels = ds.sizes['nVertLevels']

    vert_index = xr.DataArray.from_dict(
        {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)})
    mask_bottom = (vert_index < ds_mesh.maxLevelCell).T
    mask_bottom_edge = 0.5 * (mask_bottom.isel(nCells=cell0) +
                              mask_bottom.isel(nCells=cell1))

    if 'zMid' not in ds.keys():
        z_mid = compute_zmid(ds_mesh.bottomDepth, ds_mesh.maxLevelCell,
                             ds_mesh.layerThickness)
    else:
        z_mid = ds.zMid
    z_mid_edge = 0.5 * (z_mid.isel(nCells=cell0) +
                        z_mid.isel(nCells=cell1))

    mask = np.logical_and(np.logical_and(z_mid_edge >= -max_depth,
                                         z_mid_edge <= -min_depth),
                          mask_bottom_edge)
    normal_velocity = normal_velocity.where(mask)
    layer_thickness_edge = layer_thickness_edge.where(mask)
    transport = ds_mesh.dvEdge[inner_edges] * \
        (layer_thickness_edge * normal_velocity).sum(dim='nVertLevels')

    return inner_edges, transport
