import numpy as np
import xarray as xr


def compute_vertically_integrated_vorticity(
    ds_mesh,
    vert_integ_velocity,
    logger=None,
):
    """
    Compute the vertically integrated vorticity on MPAS vertices.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        A dataset containing MPAS mesh variables

    vert_integ_velocity : xarray.DataArray
        A data array containing the vertically integrated normal velocity on
        the mesh edges.

    logger : logging.Logger, optional
        A logger for the output if not stdout

    Returns
    -------
    vert_integ_vorticity : xarray.DataArray
        The vertically integrated velocity on the mesh edges

    edge_sign_on_vertex : xarray.DataArray
        The sign of the edges on the vertices
    """
    var_list = [
        'areaTriangle',
        'dcEdge',
        'edgesOnVertex',
        'verticesOnEdge',
        'latVertex',
    ]
    ds_mesh = ds_mesh[var_list].as_numpy()
    vert_integ_velocity = vert_integ_velocity.as_numpy()

    if logger:
        logger.info('  Computing edge signs on vertices.')
    # Compute the sign of edges on vertices to determine flow direction
    edge_sign_on_vertex = _compute_edge_sign_on_vertex(ds_mesh)

    if logger:
        logger.info('  Computing vertically integrated vorticity.')
    # Compute the vorticity on vertices from the integrated velocity
    vert_integ_vorticity = _compute_vert_integ_vorticity(
        ds_mesh, vert_integ_velocity, edge_sign_on_vertex
    )

    return vert_integ_vorticity, edge_sign_on_vertex


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
        valid_edge = np.logical_and(
            eov >= 0, np.logical_and(v0_on_edge >= 0, v1_on_edge >= 0)
        )

        mask = np.logical_and(valid_edge, v0_on_edge == vertices)
        edge_sign_on_vertex[mask, iedge] = -1

        mask = np.logical_and(valid_edge, v1_on_edge == vertices)
        edge_sign_on_vertex[mask, iedge] = 1

    edge_sign_on_vertex = xr.DataArray(
        edge_sign_on_vertex,
        dims=('nVertices', 'vertexDegree'),
    )
    return edge_sign_on_vertex


def _compute_vert_integ_vorticity(
    ds_mesh, vert_integ_velocity, edge_sign_on_vertex
):
    area_vertex = ds_mesh.areaTriangle
    dc_edge = ds_mesh.dcEdge
    edges_on_vertex = ds_mesh.edgesOnVertex - 1

    vertex_degree = ds_mesh.sizes['vertexDegree']

    vert_integ_vorticity = xr.zeros_like(ds_mesh.latVertex)
    for iedge in range(vertex_degree):
        eov = edges_on_vertex.isel(vertexDegree=iedge)
        edge_sign = edge_sign_on_vertex.isel(vertexDegree=iedge)
        dc = dc_edge.isel(nEdges=eov)
        vert_integ_vel = vert_integ_velocity.isel(nEdges=eov)
        vert_integ_vorticity += dc / area_vertex * edge_sign * vert_integ_vel

    return vert_integ_vorticity
