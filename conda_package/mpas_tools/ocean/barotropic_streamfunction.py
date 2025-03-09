import logging
import sys

import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import xarray as xr

from mpas_tools.ocean.depth import compute_zmid


def compute_barotropic_streamfunction(ds_mesh, ds, logger=None,
                                      min_depth=None, max_depth=None,
                                      prefix='timeMonthly_avg_',
                                      time_index=None, include_bolus=False,
                                      include_submesoscale=False):
    """
    Compute barotropic streamfunction

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        A dataset containing MPAS mesh variables

    ds : xarray.Dataset
        A dataset containing MPAS output variables ``normalVelocity`` and
        ``layerThickness`` (possibly with a ``prefix``)

    logger : logging.Logger, optional
        A logger for the output if not stdout

    min_depth : float, optional
        The minimum depth (positive up) to compute BSF over

    max_depth : float, optional
        The maximum depth (positive up) to compute BSF over

    prefix : str, optional
        The prefix on the ``normalVelocity`` and ``layerThickness`` variables

    time_index : int, optional
        The time at which to index ``ds`` (if it has ``Time`` as a dimension)

    include_bolus : bool, optional
        Whether to include the GM bolus velocity in the computation

    include_submesoscale : bool, optional
        Whether to include the submesoscale velocity in the computation

    Returns
    -------
    bsf_vertex : xarray.DataArray
        The barotropic streamfunction in Sv on vertices
    """

    useStdout = logger is None
    if useStdout:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.setLevel(logging.INFO)

    if time_index is not None:
        ds = ds.isel(Time=time_index)

    bsf_vertex = _compute_barotropic_streamfunction_vertex(
        ds_mesh, ds, prefix, include_bolus, include_submesoscale, min_depth,
        max_depth)

    return bsf_vertex


def shift_barotropic_streamfunction(bsf_vertex, lat_range, cells_on_vertex,
                                    lat_vertex):
    """
    Shift the barotropic streamfunction to be zero on average at the boundary
    over the given latitude range

    Parameters
    ----------
    bsf_vertex : xarray.DataArray
        The barotropic streamfunction in Sv on vertices

    lat_range : list of float
        The latitude range over which to set the mean boundary value of the BSF
        to zero

    cells_on_vertex : xarray.DataArray
        The zero-based cell indices on each vertex

    lat_vertex : xarray.DataArray
        The latitude of each vertex

    Returns
    -------
    bsf_shifted : xarray.DataArray
        The shifted barotropic streamfunction in Sv on vertices
    """
    is_boundary_cov = cells_on_vertex == -1
    boundary_vertices = is_boundary_cov.sum(dim='vertexDegree') > 0

    boundary_vertices = np.logical_and(
        boundary_vertices,
        lat_vertex >= np.deg2rad(lat_range[0])
    )
    boundary_vertices = np.logical_and(
        boundary_vertices,
        lat_vertex <= np.deg2rad(lat_range[1])
    )

    # convert from boolean mask to indices
    boundary_vertices = np.flatnonzero(boundary_vertices.values)

    mean_boundary_bsf = bsf_vertex.isel(nVertices=boundary_vertices).mean()

    bsf_shifted = bsf_vertex - mean_boundary_bsf

    return bsf_shifted


def _compute_vert_integ_velocity(ds_mesh, ds, prefix, include_bolus,
                                 include_submesoscale, min_depth, max_depth):

    cells_on_edge = ds_mesh.cellsOnEdge - 1
    inner_edges = np.logical_and(cells_on_edge.isel(TWO=0) >= 0,
                                 cells_on_edge.isel(TWO=1) >= 0)

    # convert from boolean mask to indices
    inner_edges = np.flatnonzero(inner_edges.values)

    cell0 = cells_on_edge.isel(nEdges=inner_edges, TWO=0)
    cell1 = cells_on_edge.isel(nEdges=inner_edges, TWO=1)
    n_vert_levels = ds.sizes['nVertLevels']

    layer_thickness = ds[f'{prefix}layerThickness']
    max_level_cell = ds_mesh.maxLevelCell - 1

    vert_index = xr.DataArray.from_dict(
        {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)})
    z_mid = compute_zmid(ds_mesh.bottomDepth, ds_mesh.maxLevelCell,
                         layer_thickness)
    z_mid_edge = 0.5*(z_mid.isel(nCells=cell0) +
                      z_mid.isel(nCells=cell1))

    normal_velocity = ds[f'{prefix}normalVelocity']
    if include_bolus:
        normal_velocity += ds[f'{prefix}normalGMBolusVelocity']
    if include_submesoscale:
        normal_velocity += ds[f'{prefix}normalMLEvelocity']
    normal_velocity = normal_velocity.isel(nEdges=inner_edges)

    layer_thickness_edge = 0.5*(layer_thickness.isel(nCells=cell0) +
                                layer_thickness.isel(nCells=cell1))
    mask_bottom = (vert_index <= max_level_cell).T
    mask_bottom_edge = np.logical_and(mask_bottom.isel(nCells=cell0),
                                      mask_bottom.isel(nCells=cell1))
    masks = [mask_bottom_edge]
    if min_depth is not None:
        masks.append(z_mid_edge <= min_depth)
    if max_depth is not None:
        masks.append(z_mid_edge >= max_depth)
    for mask in masks:
        normal_velocity = normal_velocity.where(mask)
        layer_thickness_edge = layer_thickness_edge.where(mask)

    vert_integ_velocity = np.zeros(ds_mesh.sizes['nEdges'], dtype=float)
    inner_vert_integ_vel = (
        (layer_thickness_edge * normal_velocity).sum(dim='nVertLevels'))
    vert_integ_velocity[inner_edges] = inner_vert_integ_vel.values

    vert_integ_velocity = xr.DataArray(vert_integ_velocity,
                                       dims=('nEdges',))

    return vert_integ_velocity


def _compute_edge_sign_on_vertex(ds_mesh):
    edges_on_vertex = ds_mesh.edgesOnVertex - 1
    vertices_on_edge = ds_mesh.verticesOnEdge - 1

    nvertices = ds_mesh.sizes['nVertices']
    vertex_degree = ds_mesh.sizes['vertexDegree']

    edge_sign_on_vertex = np.zeros((nvertices, vertex_degree), dtype=int)
    vertices = np.arange(nvertices)
    for iedge in range(vertex_degree):
        eov = edges_on_vertex.isel(vertexDegree=iedge)
        valid_edge = eov >= 0

        v0_on_edge = vertices_on_edge.isel(nEdges=eov, TWO=0)
        v1_on_edge = vertices_on_edge.isel(nEdges=eov, TWO=1)
        valid_edge = np.logical_and(valid_edge, v0_on_edge >= 0)
        valid_edge = np.logical_and(valid_edge, v1_on_edge >= 0)

        mask = np.logical_and(valid_edge, v0_on_edge == vertices)
        edge_sign_on_vertex[mask, iedge] = -1

        mask = np.logical_and(valid_edge, v1_on_edge == vertices)
        edge_sign_on_vertex[mask, iedge] = 1

    return edge_sign_on_vertex


def _compute_vert_integ_vorticity(ds_mesh, vert_integ_velocity,
                                  edge_sign_on_vertex):

    area_vertex = ds_mesh.areaTriangle
    dc_edge = ds_mesh.dcEdge
    edges_on_vertex = ds_mesh.edgesOnVertex - 1

    vertex_degree = ds_mesh.sizes['vertexDegree']

    vert_integ_vorticity = xr.zeros_like(ds_mesh.latVertex)
    for iedge in range(vertex_degree):
        eov = edges_on_vertex.isel(vertexDegree=iedge)
        edge_sign = edge_sign_on_vertex[:, iedge]
        dc = dc_edge.isel(nEdges=eov)
        vert_integ_vel = vert_integ_velocity.isel(nEdges=eov)
        vert_integ_vorticity += (
            dc / area_vertex * edge_sign * vert_integ_vel)

    return vert_integ_vorticity


def _compute_barotropic_streamfunction_vertex(ds_mesh, ds, prefix,
                                              include_bolus,
                                              include_submesoscale, min_depth,
                                              max_depth):
    edge_sign_on_vertex = _compute_edge_sign_on_vertex(ds_mesh)
    vert_integ_velocity = _compute_vert_integ_velocity(ds_mesh, ds, prefix,
                                                       include_bolus,
                                                       include_submesoscale,
                                                       min_depth, max_depth)
    vert_integ_vorticity = _compute_vert_integ_vorticity(
        ds_mesh, vert_integ_velocity, edge_sign_on_vertex)

    nvertices = ds_mesh.sizes['nVertices']
    vertex_degree = ds_mesh.sizes['vertexDegree']

    edges_on_vertex = ds_mesh.edgesOnVertex - 1
    vertices_on_edge = ds_mesh.verticesOnEdge - 1
    area_vertex = ds_mesh.areaTriangle
    dc_edge = ds_mesh.dcEdge
    dv_edge = ds_mesh.dvEdge

    # one equation involving vertex degree + 1 vertices for each vertex
    # plus 2 entries for the boundary condition and Lagrange multiplier
    ndata = (vertex_degree + 1) * nvertices + 2
    indices = np.zeros((2, ndata), dtype=int)
    data = np.zeros(ndata, dtype=float)

    # the laplacian on the dual mesh of the streamfunction is the
    # vertically integrated vorticity
    vertices = np.arange(nvertices, dtype=int)
    idata = (vertex_degree + 1) * vertices + 1
    indices[0, idata] = vertices
    indices[1, idata] = vertices
    for iedge in range(vertex_degree):
        eov = edges_on_vertex.isel(vertexDegree=iedge)
        dc = dc_edge.isel(nEdges=eov)
        dv = dv_edge.isel(nEdges=eov)

        v0 = vertices_on_edge.isel(nEdges=eov, TWO=0)
        v1 = vertices_on_edge.isel(nEdges=eov, TWO=1)

        edge_sign = edge_sign_on_vertex[:, iedge]

        mask = v0 == vertices
        # the difference is v1 - v0, so we want to subtract this vertex
        # when it is v0 and add it when it is v1
        this_vert_sign = np.where(mask, -1., 1.)
        # the other vertex is obviously whichever one this is not
        other_vert_index = np.where(mask, v1, v0)
        # if there are invalid vertices, we need to make sure we don't
        # index out of bounds.  The edge_sign will mask these out
        other_vert_index = np.where(other_vert_index >= 0,
                                    other_vert_index, 0)

        idata_other = idata + iedge + 1

        indices[0, idata] = vertices
        indices[1, idata] = vertices
        indices[0, idata_other] = vertices
        indices[1, idata_other] = other_vert_index

        this_data = this_vert_sign * edge_sign * dc / (dv * area_vertex)
        data[idata] += this_data
        data[idata_other] = -this_data

    # Now, the boundary condition: To begin with, we set the BSF at the
    # frist vertext to zero
    indices[0, -2] = nvertices
    indices[1, -2] = 0
    data[-2] = 1.

    # The same in the final column
    indices[0, -1] = 0
    indices[1, -1] = nvertices
    data[-1] = 1.

    # one extra spot for the Lagrange multiplier
    rhs = np.zeros(nvertices + 1, dtype=float)

    rhs[0:-1] = vert_integ_vorticity.values

    matrix = scipy.sparse.csr_matrix(
        (data, indices),
        shape=(nvertices + 1, nvertices + 1))

    solution = scipy.sparse.linalg.spsolve(matrix, rhs)

    # drop the Lagrange multiplier and convert to Sv with the desired sign
    # convention
    bsf_vertex = xr.DataArray(-1e-6 * solution[0:-1],
                              dims=('nVertices',))

    return bsf_vertex
