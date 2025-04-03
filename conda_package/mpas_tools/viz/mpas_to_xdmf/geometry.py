import numpy as np
import xarray as xr


def _build_cell_geometry(ds_mesh):
    """
    Build the geometry for fields on cells.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        The mesh dataset.

    Returns
    -------
    ds_cell_geom : xarray.Dataset
        A new Dataset containing the geometry for fields on cells.
    """

    points = [ds_mesh.xVertex, ds_mesh.yVertex, ds_mesh.zVertex]

    nverts_on_cell = ds_mesh.nEdgesOnCell
    # convert to python zero-based indexing
    vertices_on_cell = ds_mesh.verticesOnCell - 1

    max_vertices = ds_mesh.sizes['maxEdges']

    cells = []
    last = vertices_on_cell.isel(maxEdges=0)
    # fill in invalid vertices with the last in the polygon
    for ivert in range(max_vertices):
        verts = vertices_on_cell.isel(maxEdges=ivert)
        mask = ivert < nverts_on_cell
        connect = xr.where(mask, verts, last)
        last = connect
        cells.append(connect)

    ds_cell_geom = xr.Dataset()
    ds_cell_geom['points'] = xr.concat(points, dim='R3').transpose(
        'nVertices', 'R3'
    )
    ds_cell_geom['cells'] = xr.concat(cells, dim='nv').transpose(
        'nCells', 'nv'
    )
    return ds_cell_geom


def _build_edge_geometry(ds_mesh):
    """
    Build the geometry for fields on edges.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        The mesh dataset.

    Returns
    -------
    ds_edge_geom : xarray.Dataset
        A new Dataset containing the geometry for fields on edges.
    """

    points = []
    # the vertices of edge polygons are made up of both vertex and cell points
    for dim in ['x', 'y', 'z']:
        points.append(
            xr.DataArray(
                np.append(
                    ds_mesh[f'{dim}Vertex'].values,
                    ds_mesh[f'{dim}Cell'].values,
                ),
                dims=['nPoints'],
            )
        )

    nvertices = ds_mesh.sizes['nVertices']

    # convert to python zero-based indexing
    vertices_on_edge = ds_mesh.verticesOnEdge - 1
    cells_on_edge = ds_mesh.cellsOnEdge - 1

    cells = []
    v0 = vertices_on_edge.isel(TWO=0)
    v1 = vertices_on_edge.isel(TWO=1)
    c0 = cells_on_edge.isel(TWO=0)
    c1 = cells_on_edge.isel(TWO=1)

    # if cells are valid, offset them by nvertices because that's how we
    # concatinated the vertices together above; if invalid, reuse the previous
    # vertex
    c0 = xr.where(c0 >= 0, c0 + nvertices, v1)
    c1 = xr.where(c1 >= 0, c1 + nvertices, v0)
    cells = [c0, v0, c1, v1]

    ds_edge_geom = xr.Dataset()
    ds_edge_geom['points'] = xr.concat(points, dim='R3').transpose(
        'nPoints', 'R3'
    )
    ds_edge_geom['cells'] = xr.concat(cells, dim='nv').transpose(
        'nEdges', 'nv'
    )
    return ds_edge_geom


def _build_vertex_geometry(ds_mesh):
    """
    Build the geometry for fields on vertices.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        The mesh dataset.

    Returns
    -------
    ds_vertex_geom : xarray.Dataset
        A new Dataset containing the geometry for fields on vertices.
    """

    points = []
    # the vertices of vertex polygons potentially involve vertex, edge, and
    # cell points
    for dim in ['x', 'y', 'z']:
        points.append(
            xr.DataArray(
                np.concatenate(
                    [
                        ds_mesh[f'{dim}Cell'].values,
                        ds_mesh[f'{dim}Edge'].values,
                        ds_mesh[f'{dim}Vertex'].values,
                    ]
                ),
                dims=['nPoints'],
            )
        )

    nvertices = ds_mesh.sizes['nVertices']
    nedges = ds_mesh.sizes['nEdges']
    ncells = ds_mesh.sizes['nCells']

    # convert to python zero-based indexing
    cells_on_vertex = ds_mesh.cellsOnVertex - 1
    edges_on_vertex = ds_mesh.edgesOnVertex - 1

    vertex_degree = ds_mesh.sizes['vertexDegree']

    cells = []
    # offset because of how we concatenated the vertices above
    verts = xr.DataArray(
        np.arange(nvertices) + ncells + nedges, dims=['nVertices']
    )
    # fill in invalid cells with the last in the polygon
    for ivert in range(vertex_degree):
        eov = edges_on_vertex.isel(vertexDegree=ivert)
        mask = eov >= 0
        # offset the edge index by ncells because of how we concatenated the
        # vertices above
        cells.append(xr.where(mask, eov + ncells, verts))
        cov = cells_on_vertex.isel(vertexDegree=ivert)
        mask = cov >= 0
        cells.append(xr.where(mask, cov, verts))

    ds_vertex_geom = xr.Dataset()
    ds_vertex_geom['points'] = xr.concat(points, dim='R3').transpose(
        'nPoints', 'R3'
    )
    ds_vertex_geom['cells'] = xr.concat(cells, dim='nv').transpose(
        'nVertices', 'nv'
    )
    return ds_vertex_geom
