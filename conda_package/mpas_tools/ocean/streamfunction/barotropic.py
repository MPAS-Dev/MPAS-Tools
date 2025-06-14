import logging
import os
import sys

import networkx as nx
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import xarray as xr

from mpas_tools.io import write_netcdf
from mpas_tools.ocean.streamfunction.velocity import (
    compute_vertically_integrated_velocity,
)
from mpas_tools.ocean.streamfunction.vorticity import (
    compute_vertically_integrated_vorticity,
)


def compute_barotropic_streamfunction(
    ds_mesh,
    ds,
    logger=None,
    min_depth=None,
    max_depth=None,
    prefix='timeMonthly_avg_',
    time_index=None,
    include_bolus=False,
    include_submesoscale=False,
    quiet=False,
    horiz_chunk=10000,
    tmp_dir=None,
):
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

    quiet : bool, optional
        If True, suppress all logging output
        If False, log all output to the logger

    horiz_chunk : int, optional
        The number of edges to chunk the dataset by when computing
        the vertically integrated velocity. This is useful for
        large datasets to avoid memory issues. Set this to ``None`` to disable
        chunking.

    tmp_dir : str, optional
        A temporary directory to use for intermediate files. This is useful
        for large datasets to avoid memory issues. If None, no temporary
        directory is used.

    Returns
    -------
    bsf_vertex : xarray.DataArray
        The barotropic streamfunction in Sv on vertices
    """

    if quiet:
        logger = None
    elif logger is None:
        logger = logging.getLogger()
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.setLevel(logging.INFO)

    if time_index is None:
        if 'Time' in ds.dims:
            raise ValueError(
                'time_index must be provided if "Time" is a dimension of ds'
            )
    else:
        ds = ds.isel(Time=time_index)

    if logger:
        logger.info('Computing barotropic streamfunction.')

    bsf_vertex = _compute_barotropic_streamfunction_vertex(
        ds_mesh,
        ds,
        prefix,
        include_bolus,
        include_submesoscale,
        min_depth,
        max_depth,
        logger,
        horiz_chunk,
        tmp_dir,
    )

    if logger:
        logger.info('  Done.')

    return bsf_vertex


def shift_barotropic_streamfunction(
    bsf_vertex, lat_range, cells_on_vertex, lat_vertex, logger=None
):
    """
    Shift the barotropic streamfunction to be zero on average at the boundary
    over the given latitude range

    Parameters
    ----------
    bsf_vertex : xarray.DataArray
        The barotropic streamfunction in Sv on vertices

    lat_range : list of float
        The latitude range in degrees over which to set the mean boundary value
        of the BSF to zero

    cells_on_vertex : xarray.DataArray
        The zero-based cell indices on each vertex

    lat_vertex : xarray.DataArray
        The latitude of each vertex in radians

    logger : logging.Logger, optional
        A logger for the output

    Returns
    -------
    bsf_shifted : xarray.DataArray
        The shifted barotropic streamfunction in Sv on vertices
    """
    is_boundary_cov = cells_on_vertex == -1
    boundary_vertices = np.logical_and(
        is_boundary_cov.sum(dim='vertexDegree') > 0,
        np.logical_and(
            lat_vertex >= np.deg2rad(lat_range[0]),
            lat_vertex <= np.deg2rad(lat_range[1]),
        ),
    )

    # convert from boolean mask to indices
    boundary_vertices = np.flatnonzero(boundary_vertices.values)

    mean_boundary_bsf = bsf_vertex.isel(nVertices=boundary_vertices).mean()

    if logger:
        logger.info(
            f'    Mean BSF on boundary vertices in range {lat_range} '
            f'is {mean_boundary_bsf.values:.4f} Sv'
        )

    bsf_shifted = bsf_vertex - mean_boundary_bsf

    return bsf_shifted


def _build_minimal_boundary_constraints(
    boundary_vertex0, boundary_vertex1, logger
):
    """
    Construct a minimal set of boundary constraints that tie each connected
    loop together without introducing redundancy.

    Parameters
    ----------
    boundary_vertex0 : xarray.DataArray
        The first vertex in each unique pair of boundary vertices.

    boundary_vertex1 : xarray.DataArray
        The second vertex in each unique pair of boundary vertices.

    logger : logging.Logger, optional
        Logger for logging messages.

    Returns
    -------
    minimal_constraints : list of tuple
        A minimal set of edges (vertex pairs) that constrain the boundary.
    """
    # Create a graph from the boundary edges
    graph = nx.Graph()
    edges = list(zip(boundary_vertex0.values, boundary_vertex1.values))
    graph.add_edges_from(edges)

    minimal_constraints = []

    # Loop over connected components (disjoint loops)
    connected_components = list(nx.connected_components(graph))
    if logger:
        logger.info(
            f'    Found {len(connected_components)} independent boundary '
            f'loops.'
        )

    for component in connected_components:
        subgraph = graph.subgraph(component)

        # Find a spanning tree for the component
        spanning_tree = nx.minimum_spanning_tree(subgraph)

        # Add the edges of the spanning tree to the constraints
        minimal_constraints.extend(spanning_tree.edges)

    return minimal_constraints


def _identify_boundary_vertices(ds_mesh, logger, all_vertices):
    """
    Identify boundary vertices and edges in the mesh, and remove redundant
    vertices that lead to overdetermined constraints.
    """

    if logger:
        logger.info('  Identifying boundary vertices.')
    # Identify boundary vertices
    nedges = ds_mesh.sizes['nEdges']
    cells_on_vertex = ds_mesh.cellsOnVertex - 1
    cells_on_edge = ds_mesh.cellsOnEdge - 1
    vertices_on_edge = ds_mesh.verticesOnEdge - 1

    boundary_mask = (cells_on_vertex == -1).any(dim='vertexDegree')

    boundary_vertices = all_vertices.where(boundary_mask, drop=True).astype(
        int
    )

    boundary_edge_mask = (cells_on_edge == -1).any(dim='TWO')
    all_edges = xr.DataArray(np.arange(nedges, dtype=int), dims=('nEdges',))
    boundary_edges = all_edges.where(boundary_edge_mask, drop=True).astype(int)

    boundary_vertex0 = vertices_on_edge.isel(nEdges=boundary_edges, TWO=0)
    boundary_vertex1 = vertices_on_edge.isel(nEdges=boundary_edges, TWO=1)

    if logger:
        logger.info('  Detect boundary loops and remove redundant pairs.')
    # find each independent loop of boundary vertices and determine a set of
    # pairs that does not close the loop (i.e. a spanning tree), avoiding
    # overdetermined constraints
    minimal_constraints = _build_minimal_boundary_constraints(
        boundary_vertex0, boundary_vertex1, logger
    )

    # Unpack minimal_constraints into boundary_vertex0 and boundary_vertex1
    boundary_vertex0, boundary_vertex1 = zip(*minimal_constraints)
    boundary_vertex0 = xr.DataArray(
        np.array(boundary_vertex0), dims=('nVertices',)
    )
    boundary_vertex1 = xr.DataArray(
        np.array(boundary_vertex1), dims=('nVertices',)
    )

    nboundary = boundary_vertex0.sizes['nVertices']

    if logger:
        nboundary_removed = boundary_vertices.sizes['nVertices'] - nboundary
        logger.info(
            f'  Removed {nboundary_removed} redundant boundary vertices.'
        )

    return boundary_vertex0, boundary_vertex1, nboundary


def _assemble_matrix(
    ds_mesh, logger, edge_sign_on_vertex, vert_integ_vorticity
):
    """
    Assemble the sparse matrix for the Poisson equation.
    """

    var_list = [
        'edgesOnVertex',
        'verticesOnEdge',
        'areaTriangle',
        'dcEdge',
        'dvEdge',
        'cellsOnVertex',
        'cellsOnEdge',
    ]
    ds_mesh = ds_mesh[var_list].as_numpy()
    edge_sign_on_vertex = edge_sign_on_vertex.as_numpy()
    vert_integ_vorticity = vert_integ_vorticity.as_numpy()

    nvertices = ds_mesh.sizes['nVertices']
    vertex_degree = ds_mesh.sizes['vertexDegree']
    edges_on_vertex = ds_mesh.edgesOnVertex - 1
    vertices_on_edge = ds_mesh.verticesOnEdge - 1
    area_vertex = ds_mesh.areaTriangle
    dc_edge = ds_mesh.dcEdge
    dv_edge = ds_mesh.dvEdge

    all_vertices = xr.DataArray(
        np.arange(nvertices, dtype=int), dims=('nVertices',)
    )

    boundary_vertex0, boundary_vertex1, nboundary = (
        _identify_boundary_vertices(ds_mesh, logger, all_vertices)
    )

    if logger:
        logger.info('  Assembling sparse matrix for the Poisson equation.')
    # Assemble the sparse matrix for solving the Poisson equation:
    #   * the Poisson equation at each vertex involves vertex degree + 1 terms
    #   * the boundary conditions involve 2 vertices and are duplicated in
    #     the form of a Lagrange multiplier
    #   * the unique solution is ensured by adding a constraint on the
    #     streamfunction at the first vertex, again as a Lagrange multiplier
    ndata = (vertex_degree + 1) * nvertices + 4 * nboundary + 2
    nmatrix = nvertices + nboundary + 1
    indices = np.zeros((2, ndata), dtype=int)
    data = np.zeros(ndata, dtype=float)

    # Fill the Poisson equation for the BSF at each vertex will be equal to the
    # vertically integrated vorticity
    idata = (vertex_degree + 1) * all_vertices.values
    rows = all_vertices.values
    indices[0, idata] = rows
    indices[1, idata] = all_vertices.values
    for iedge in range(vertex_degree):
        eov = edges_on_vertex.isel(vertexDegree=iedge)
        dc = dc_edge.isel(nEdges=eov)
        dv = dv_edge.isel(nEdges=eov)

        v0 = vertices_on_edge.isel(nEdges=eov, TWO=0)
        v1 = vertices_on_edge.isel(nEdges=eov, TWO=1)

        edge_sign = edge_sign_on_vertex.isel(vertexDegree=iedge)

        mask = v0 == all_vertices
        this_vert_sign = xr.where(mask, -1.0, 1.0)
        other_vert_index = xr.where(mask, v1, v0)
        other_vert_index = xr.where(other_vert_index >= 0, other_vert_index, 0)

        idata_other = idata + iedge + 1

        indices[0, idata_other] = rows
        indices[1, idata_other] = other_vert_index.values

        this_data = this_vert_sign * edge_sign * dc / (dv * area_vertex)
        data[idata] += this_data.values
        data[idata_other] = -this_data.values

    # Add boundary conditions to the matrix
    #  The difference in the BSF between adjacent boundary vertices is
    #  zero
    idata = (vertex_degree + 1) * nvertices + 2 * np.arange(nboundary)
    rows = nvertices + np.arange(nboundary)
    indices[0, idata] = rows
    indices[1, idata] = boundary_vertex0.values
    data[idata] = -1.0

    idata += 1
    indices[0, idata] = rows
    indices[1, idata] = boundary_vertex1.values
    data[idata] = 1.0

    # Now the transpose
    idata = (
        (vertex_degree + 1) * nvertices
        + 2 * nboundary
        + 2 * np.arange(nboundary)
    )
    col = nvertices + np.arange(nboundary)
    indices[0, idata] = boundary_vertex0.values
    indices[1, idata] = col
    data[idata] = -1.0

    idata += 1
    indices[0, idata] = boundary_vertex1.values
    indices[1, idata] = col
    data[idata] = 1.0

    # Add gauge constraints to ensure a unique solution
    # The BSF at vertex 0 will be zero
    idata = ndata - 2
    indices[0, idata] = nmatrix - 1
    indices[1, idata] = 0
    data[idata] = 1.0

    # And the transpose
    idata = ndata - 1
    indices[0, idata] = 0
    indices[1, idata] = nmatrix - 1
    data[idata] = 1.0

    # Assemble the right-hand side of the equation:
    #  * the vertically integrated vorticity at each vertex
    #  * the boundary conditions are zero
    #  * the gauge constraints are zero
    rhs = np.zeros(nmatrix, dtype=float)
    rhs[0:nvertices] = vert_integ_vorticity.values

    return indices, data, rhs, nmatrix


def _compute_barotropic_streamfunction_vertex(
    ds_mesh,
    ds,
    prefix,
    include_bolus,
    include_submesoscale,
    min_depth,
    max_depth,
    logger,
    horiz_chunk,
    tmp_dir,
):
    """
    Compute the barotropic streamfunction on vertices.

    This function solves a Poisson equation to compute the barotropic
    streamfunction, which integrates vertically integrated velocity
    divergence to obtain the streamfunction.
    """
    if logger:
        logger.info('  Computing vertically integrated velocity.')
    # Compute the vertically integrated velocity on edges
    vert_integ_velocity = compute_vertically_integrated_velocity(
        ds_mesh=ds_mesh,
        ds=ds,
        logger=logger,
        min_depth=min_depth,
        max_depth=max_depth,
        prefix=prefix,
        include_bolus=include_bolus,
        include_submesoscale=include_submesoscale,
        nedges_chunk=horiz_chunk,
    )

    if tmp_dir is not None:
        if logger:
            logger.info(
                '  Writing out and reading in vertically integrated velocity.'
            )
        # write out and read back the vertically integrated velocity
        ds_out = xr.Dataset()
        ds_out['vertIntegNormalVelocity'] = vert_integ_velocity
        filename = os.path.join(tmp_dir, 'bsf_vert_integ_normal_vel.nc')
        write_netcdf(ds_out, filename, logger=logger)
        ds_in = xr.open_dataset(filename)
        vert_integ_velocity = ds_in['vertIntegNormalVelocity']

    vert_integ_vorticity, edge_sign_on_vertex = (
        compute_vertically_integrated_vorticity(
            ds_mesh=ds_mesh,
            vert_integ_velocity=vert_integ_velocity,
            logger=logger,
        )
    )

    if tmp_dir is not None:
        if logger:
            logger.info(
                '  Writing out and reading in vertically integrated vorticity.'
            )
        # write out and read back the vertically integrated velocity
        ds_out = xr.Dataset()
        ds_out['vertIntegVorticity'] = vert_integ_vorticity
        ds_out['edgeSignOnVertex'] = edge_sign_on_vertex
        filename = os.path.join(tmp_dir, 'bsf_vert_integ_vorticity.nc')
        write_netcdf(ds_out, filename, logger=logger)
        ds_in = xr.open_dataset(filename)
        vert_integ_vorticity = ds_in['vertIntegVorticity']
        edge_sign_on_vertex = ds_in['edgeSignOnVertex']

    indices, data, rhs, nmatrix = _assemble_matrix(
        ds_mesh,
        logger,
        edge_sign_on_vertex,
        vert_integ_vorticity,
    )

    if logger:
        logger.info('  Solving the sparse linear system.')
    # Solve the sparse linear system
    matrix = scipy.sparse.csr_matrix((data, indices), shape=(nmatrix, nmatrix))
    solution = scipy.sparse.linalg.spsolve(matrix, rhs)

    if logger:
        logger.info('  Finalizing the barotropic streamfunction.')
    # Convert the solution to the barotropic streamfunction
    nvertices = ds_mesh.sizes['nVertices']
    bsf_vertex = xr.DataArray(
        -1e-6 * solution[0:nvertices], dims=('nVertices',)
    )

    return bsf_vertex
