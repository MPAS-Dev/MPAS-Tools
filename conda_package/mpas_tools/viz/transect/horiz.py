import numpy as np
import xarray as xr
from scipy.spatial import KDTree
from shapely.geometry import LineString, Point

from mpas_tools.transects import (
    cartesian_to_lon_lat,
    lon_lat_to_cartesian,
    subdivide_great_circle,
    subdivide_planar,
)
from mpas_tools.vector import Vector


def mesh_to_triangles(ds_mesh):
    """
    Construct a dataset in which each MPAS cell is divided into the triangles
    connecting pairs of adjacent vertices to cell centers.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        An MPAS mesh

    Returns
    -------
    ds_tris : xarray.Dataset
        A dataset that defines triangles connecting pairs of adjacent vertices
        to cell centers as well as the cell index that each triangle is in and
        cell indices and weights for interpolating data defined at cell centers
        to triangle nodes.  ``ds_tris`` includes variables ``triCellIndices``,
        the cell that each triangle is part of; ``nodeCellIndices`` and
        ``nodeCellWeights``, the indices and weights used to interpolate from
        MPAS cell centers to triangle nodes; Cartesian coordinates ``xNode``,
        ``yNode``, and ``zNode``; and ``lonNode``` and ``latNode`` in radians.
        ``lonNode`` is guaranteed to be within 180 degrees of the cell center
        corresponding to ``triCellIndices``.  Nodes always have a
        counterclockwise winding.

    """
    n_vertices_on_cell = ds_mesh.nEdgesOnCell.values
    vertices_on_cell = ds_mesh.verticesOnCell.values - 1
    cells_on_vertex = ds_mesh.cellsOnVertex.values - 1

    on_a_sphere = ds_mesh.attrs['on_a_sphere'].strip() == 'YES'
    is_periodic = False
    x_period = None
    y_period = None
    if not on_a_sphere:
        is_periodic = ds_mesh.attrs['is_periodic'].strip() == 'YES'
        if is_periodic:
            x_period = ds_mesh.attrs['x_period']
            y_period = ds_mesh.attrs['y_period']

    kite_areas_on_vertex = ds_mesh.kiteAreasOnVertex.values

    n_triangles = np.sum(n_vertices_on_cell)

    max_edges = ds_mesh.sizes['maxEdges']
    n_cells = ds_mesh.sizes['nCells']
    if ds_mesh.sizes['vertexDegree'] != 3:
        raise ValueError(
            'mesh_to_triangles only supports meshes with vertexDegree = 3'
        )

    # find the third vertex for each triangle
    next_vertex = -1 * np.ones(vertices_on_cell.shape, int)
    for i_vertex in range(max_edges):
        valid = i_vertex < n_vertices_on_cell
        invalid = np.logical_not(valid)
        vertices_on_cell[invalid, i_vertex] = -1
        nv = n_vertices_on_cell[valid]
        cell_indices = np.arange(0, n_cells)[valid]
        i_next = np.where(i_vertex < nv - 1, i_vertex + 1, 0)
        next_vertex[:, i_vertex][valid] = vertices_on_cell[
            cell_indices, i_next
        ]

    valid = vertices_on_cell >= 0
    vertices_on_cell = vertices_on_cell[valid]
    next_vertex = next_vertex[valid]

    # find the cell index for each triangle
    tri_cell_indices, _ = np.meshgrid(
        np.arange(0, n_cells), np.arange(0, max_edges), indexing='ij'
    )
    tri_cell_indices = tri_cell_indices[valid]

    # find list of cells and weights for each triangle node
    node_cell_indices = -1 * np.ones((n_triangles, 3, 3), dtype=int)
    node_cell_weights = np.zeros((n_triangles, 3, 3))

    # the first node is at the cell center, so the value is just the one from
    # that cell
    node_cell_indices[:, 0, 0] = tri_cell_indices
    node_cell_weights[:, 0, 0] = 1.0

    # the other 2 nodes are associated with vertices
    node_cell_indices[:, 1, :] = cells_on_vertex[vertices_on_cell, :]
    node_cell_weights[:, 1, :] = kite_areas_on_vertex[vertices_on_cell, :]
    node_cell_indices[:, 2, :] = cells_on_vertex[next_vertex, :]
    node_cell_weights[:, 2, :] = kite_areas_on_vertex[next_vertex, :]

    weight_sum = np.sum(node_cell_weights, axis=2)
    for i_node in range(3):
        node_cell_weights[:, :, i_node] = (
            node_cell_weights[:, :, i_node] / weight_sum
        )

    ds_tris = xr.Dataset()
    ds_tris['triCellIndices'] = ('nTriangles', tri_cell_indices)
    ds_tris['nodeCellIndices'] = (
        ('nTriangles', 'nNodes', 'nInterp'),
        node_cell_indices,
    )
    ds_tris['nodeCellWeights'] = (
        ('nTriangles', 'nNodes', 'nInterp'),
        node_cell_weights,
    )

    # get Cartesian and lon/lat coordinates of each node
    for prefix in ['x', 'y', 'z', 'lat', 'lon']:
        out_var = f'{prefix}Node'
        cell_var = f'{prefix}Cell'
        vertex_var = f'{prefix}Vertex'
        coord = np.zeros((n_triangles, 3))
        coord[:, 0] = ds_mesh[cell_var].values[tri_cell_indices]
        coord[:, 1] = ds_mesh[vertex_var].values[vertices_on_cell]
        coord[:, 2] = ds_mesh[vertex_var].values[next_vertex]
        ds_tris[out_var] = (('nTriangles', 'nNodes'), coord)

    # nothing obvious we can do about triangles containing the poles

    if on_a_sphere:
        ds_tris = _fix_periodic_tris(
            ds_tris, periodic_var='lonNode', period=2 * np.pi
        )
    elif is_periodic:
        ds_tris = _fix_periodic_tris(
            ds_tris, periodic_var='xNode', period=x_period
        )
        ds_tris = _fix_periodic_tris(
            ds_tris, periodic_var='yNode', period=y_period
        )

    return ds_tris


def make_triangle_tree(ds_tris):
    """
    Make a KD-Tree for finding triangle edges that are near enough to transect
    segments that they might intersect

    Parameters
    ----------
    ds_tris : xarray.Dataset
        A dataset that defines triangles, the results of calling
        :py:func:`mpas_tools.viz.transect.horiz.mesh_to_triangles()`

    Returns
    -------
    tree : scipy.spatial.KDTree
        A tree of edge centers from triangles making up an MPAS mesh
    """

    n_triangles = ds_tris.sizes['nTriangles']
    n_nodes = ds_tris.sizes['nNodes']
    node_coords = np.zeros((n_triangles * n_nodes, 3))
    node_coords[:, 0] = ds_tris.xNode.values.ravel()
    node_coords[:, 1] = ds_tris.yNode.values.ravel()
    node_coords[:, 2] = ds_tris.zNode.values.ravel()

    next_tri, next_node = np.meshgrid(
        np.arange(n_triangles),
        np.mod(np.arange(n_nodes) + 1, 3),
        indexing='ij',
    )
    nextIndices = n_nodes * next_tri.ravel() + next_node.ravel()

    # edge centers are half way between adjacent nodes (ignoring great-circle
    # distance)
    edgeCoords = 0.5 * (node_coords + node_coords[nextIndices, :])

    tree = KDTree(data=edgeCoords, copy_data=True)
    return tree


def find_spherical_transect_cells_and_weights(
    lon_transect,
    lat_transect,
    ds_tris,
    ds_mesh,
    tree,
    degrees=True,
    earth_radius=None,
    subdivision_res=10e3,
):
    """
    Find "nodes" where the transect intersects the edges of the triangles
    that make up MPAS cells.

    Parameters
    ----------
    lon_transect : xarray.DataArray
        The longitude of segments making up the transect

    lat_transect : xarray.DataArray
        The latitude of segments making up the transect

    ds_tris : xarray.Dataset
        A dataset that defines triangles, the results of calling
        :py:func:`mpas_tools.viz.transect.horiz.mesh_to_triangles()`

    ds_mesh : xarray.Dataset
        A data set with the full MPAS mesh.

    tree : scipy.spatial.KDTree
        A tree of edge centers from triangles making up an MPAS mesh, the
        return value from
        :py:func:`mpas_tools.viz.transect.horiz.make_triangle_tree()`

    degrees : bool, optional
        Whether ``lon_transect`` and ``lat_transect`` are in degrees (as
        opposed to radians).

    subdivision_res : float, optional
        Resolution in m to use to subdivide the transect when looking for
        intersection candidates.  Should be small enough that curvature is
        small.

    earth_radius : float, optional
        The radius of the Earth in meters, taken from the `sphere_radius`
        global attribute if not provided

    Returns
    -------
    ds_out : xarray.Dataset
        A dataset that contains "nodes" where the transect intersects the
        edges of the triangles in ``ds_tris``.  The nodes also includes the two
        end points of the transect, which typically lie within triangles. Each
        internal node (that is, not including the end points) is purposefully
        repeated twice, once for each triangle that node touches.  This allows
        for discontinuous fields between triangles (e.g. if one wishes to plot
        constant values on each MPAS cell).  The Cartesian and lon/lat
        coordinates of these nodes are ``xCartNode``, ``yCartNode``,
        ``zCartNode``, ``lonNode`` and ``latNode``.  The distance along the
        transect of each intersection is ``dNode``. The index of the triangle
        and the first triangle node in ``ds_tris`` associated with each
        intersection node are given by ``horizTriangleIndices`` and
        ``horizTriangleNodeIndices``, respectively. The second node on the
        triangle for the edge associated with the intersection is given by
        ``numpy.mod(horizTriangleNodeIndices + 1, 3)``.

        The MPAS cell that a given node belongs to is given by
        ``horizCellIndices``. Each node also has an associated set of 6
        ``interpHorizCellIndices`` and ``interpHorizCellWeights`` that can be
        used to interpolate from MPAS cell centers to nodes first with
        area-weighted averaging to MPAS vertices and then linear interpolation
        along triangle edges.  Some of the weights may be zero, in which case
        the associated ``interpHorizCellIndices`` will be -1.

        Finally, ``lonTransect`` and ``latTransect`` are included in the
        dataset, along with Cartesian coordinates ``xCartTransect``,
        ``yCartTransect``, `zCartTransect``, and ``dTransect``, the
        great-circle distance along the transect of each original transect
        point.  In order to interpolate values (e.g. observations) from the
        original transect points to the intersection nodes, linear
        interpolation indices ``transectIndicesOnHorizNode`` and weights
        ``transectWeightsOnHorizNode`` are provided.  The values at nodes are
        found by::

          nodeValues = ((transectValues[transectIndicesOnHorizNode] *
                         transectWeightsOnHorizNode)
                        + (transectValues[transectIndicesOnHorizNode+1] *
                           (1.0 - transectWeightsOnHorizNode))
    """
    if earth_radius is None:
        earth_radius = ds_mesh.attrs['sphere_radius']
    buffer = np.maximum(
        np.amax(ds_mesh.dvEdge.values), np.amax(ds_mesh.dcEdge.values)
    )

    x, y, z = lon_lat_to_cartesian(
        lon_transect, lat_transect, earth_radius, degrees
    )

    n_nodes = ds_tris.sizes['nNodes']
    node_cell_weights = ds_tris.nodeCellWeights.values
    node_cell_indices = ds_tris.nodeCellIndices.values

    x_node = ds_tris.xNode.values.ravel()
    y_node = ds_tris.yNode.values.ravel()
    z_node = ds_tris.zNode.values.ravel()

    d_transect = np.zeros(lon_transect.shape)

    d_node = None
    x_out = None
    y_out = None
    z_out = None
    tris = None
    nodes = None
    interp_cells = None
    cell_weights = None

    n_horiz_weights = 6

    first = True

    d_start = 0.0
    for seg_index in range(len(x) - 1):
        transectv0 = Vector(
            x[seg_index].values, y[seg_index].values, z[seg_index].values
        )
        transectv1 = Vector(
            x[seg_index + 1].values,
            y[seg_index + 1].values,
            z[seg_index + 1].values,
        )

        sub_slice = slice(seg_index, seg_index + 2)
        x_sub, y_sub, z_sub, _, _ = subdivide_great_circle(
            x[sub_slice].values,
            y[sub_slice].values,
            z[sub_slice].values,
            subdivision_res,
            earth_radius,
        )

        coords = np.zeros((len(x_sub), 3))
        coords[:, 0] = x_sub
        coords[:, 1] = y_sub
        coords[:, 2] = z_sub
        radius = buffer + subdivision_res

        index_list = tree.query_ball_point(x=coords, r=radius)

        unique_indices = set()
        for indices in index_list:
            unique_indices.update(indices)

        n0_indices_cand = np.array(list(unique_indices))

        if len(n0_indices_cand) == 0:
            continue

        tris_cand = n0_indices_cand // n_nodes
        next_node_index = np.mod(n0_indices_cand + 1, n_nodes)
        n1_indices_cand = n_nodes * tris_cand + next_node_index

        n0_cand = Vector(
            x_node[n0_indices_cand],
            y_node[n0_indices_cand],
            z_node[n0_indices_cand],
        )
        n1_cand = Vector(
            x_node[n1_indices_cand],
            y_node[n1_indices_cand],
            z_node[n1_indices_cand],
        )

        intersect = Vector.intersects(n0_cand, n1_cand, transectv0, transectv1)

        n0_inter = Vector(
            n0_cand.x[intersect], n0_cand.y[intersect], n0_cand.z[intersect]
        )
        n1_inter = Vector(
            n1_cand.x[intersect], n1_cand.y[intersect], n1_cand.z[intersect]
        )

        tris_inter = tris_cand[intersect]
        n0_indices_inter = n0_indices_cand[intersect]
        n1_indices_inter = n1_indices_cand[intersect]

        intersections = Vector.intersection(
            n0_inter, n1_inter, transectv0, transectv1
        )
        intersections = Vector(
            earth_radius * intersections.x,
            earth_radius * intersections.y,
            earth_radius * intersections.z,
        )

        angular_distance = transectv0.angular_distance(intersections)

        d_node_local = d_start + earth_radius * angular_distance

        d_start += earth_radius * transectv0.angular_distance(transectv1)

        node0_inter = np.mod(n0_indices_inter, n_nodes)
        node1_inter = np.mod(n1_indices_inter, n_nodes)

        node_weights = intersections.angular_distance(
            n1_inter
        ) / n0_inter.angular_distance(n1_inter)

        weights = np.zeros((len(tris_inter), n_horiz_weights))
        cell_indices = np.zeros((len(tris_inter), n_horiz_weights), int)
        for index in range(3):
            weights[:, index] = (
                node_weights
                * node_cell_weights[tris_inter, node0_inter, index]
            )
            cell_indices[:, index] = node_cell_indices[
                tris_inter, node0_inter, index
            ]
            weights[:, index + 3] = (1.0 - node_weights) * node_cell_weights[
                tris_inter, node1_inter, index
            ]
            cell_indices[:, index + 3] = node_cell_indices[
                tris_inter, node1_inter, index
            ]

        if first:
            x_out = intersections.x
            y_out = intersections.y
            z_out = intersections.z
            d_node = d_node_local

            tris = tris_inter
            nodes = node0_inter
            interp_cells = cell_indices
            cell_weights = weights
            first = False
        else:
            x_out = np.append(x_out, intersections.x)
            y_out = np.append(y_out, intersections.y)
            z_out = np.append(z_out, intersections.z)
            d_node = np.append(d_node, d_node_local)

            tris = np.concatenate((tris, tris_inter))
            nodes = np.concatenate((nodes, node0_inter))
            interp_cells = np.concatenate((interp_cells, cell_indices), axis=0)
            cell_weights = np.concatenate((cell_weights, weights), axis=0)

        d_transect[seg_index + 1] = d_start

    epsilon = 1e-6 * subdivision_res
    (
        d_node,
        x_out,
        y_out,
        z_out,
        seg_tris,
        seg_nodes,
        interp_cells,
        cell_weights,
        valid_nodes,
    ) = _sort_intersections(
        d_node,
        tris,
        nodes,
        x_out,
        y_out,
        z_out,
        interp_cells,
        cell_weights,
        epsilon,
    )

    lon_out, lat_out = cartesian_to_lon_lat(
        x_out, y_out, z_out, earth_radius, degrees
    )

    valid_segs = seg_tris >= 0
    cell_indices = -1 * np.ones(seg_tris.shape, dtype=int)
    cell_indices[valid_segs] = ds_tris.triCellIndices.values[
        seg_tris[valid_segs]
    ]

    ds_out = xr.Dataset()
    ds_out['xCartNode'] = (('nNodes',), x_out)
    ds_out['yCartNode'] = (('nNodes',), y_out)
    ds_out['zCartNode'] = (('nNodes',), z_out)
    ds_out['dNode'] = (('nNodes',), d_node)
    ds_out['lonNode'] = (('nNodes',), lon_out)
    ds_out['latNode'] = (('nNodes',), lat_out)

    ds_out['horizTriangleIndices'] = ('nSegments', seg_tris)
    ds_out['horizCellIndices'] = ('nSegments', cell_indices)
    ds_out['horizTriangleNodeIndices'] = (
        ('nSegments', 'nHorizBounds'),
        seg_nodes,
    )
    ds_out['interpHorizCellIndices'] = (
        ('nNodes', 'nHorizWeights'),
        interp_cells,
    )
    ds_out['interpHorizCellWeights'] = (
        ('nNodes', 'nHorizWeights'),
        cell_weights,
    )
    ds_out['validNodes'] = (('nNodes',), valid_nodes)

    transect_indices_on_horiz_node = np.zeros(d_node.shape, dtype=int)
    transect_weights_on_horiz_node = np.zeros(d_node.shape)
    for trans_index in range(len(d_transect) - 1):
        d0 = d_transect[trans_index]
        d1 = d_transect[trans_index + 1]
        mask = np.logical_and(d_node >= d0, d_node < d1)
        transect_indices_on_horiz_node[mask] = trans_index
        transect_weights_on_horiz_node[mask] = (d1 - d_node[mask]) / (d1 - d0)
    # last index will get missed by the mask and needs to be handled as a
    # special case
    transect_indices_on_horiz_node[-1] = len(d_transect) - 2
    transect_weights_on_horiz_node[-1] = 0.0

    ds_out['lonTransect'] = lon_transect
    ds_out['latTransect'] = lat_transect
    ds_out['xCartTransect'] = x
    ds_out['yCartTransect'] = y
    ds_out['zCartTransect'] = z
    ds_out['dTransect'] = (lon_transect.dims, d_transect)
    ds_out['transectIndicesOnHorizNode'] = (
        ('nNodes',),
        transect_indices_on_horiz_node,
    )
    ds_out['transectWeightsOnHorizNode'] = (
        ('nNodes',),
        transect_weights_on_horiz_node,
    )

    return ds_out


def find_planar_transect_cells_and_weights(
    x_transect, y_transect, ds_tris, ds_mesh, tree, subdivision_res=10e3
):
    """
    Find "nodes" where the transect intersects the edges of the triangles
    that make up MPAS cells.

    Parameters
    ----------
    x_transect : xarray.DataArray
        The x points defining segments making up the transect

    y_transect : xarray.DataArray
        The y points defining segments making up the transect

    ds_tris : xarray.Dataset
        A dataset that defines triangles, the results of calling
        :py:func:`mpas_tools.viz.transect.horiz.mesh_to_triangles()`

    ds_mesh : xarray.Dataset
        A data set with the full MPAS mesh.

    tree : scipy.spatial.KDTree
        A tree of edge centers from triangles making up an MPAS mesh, the
        return value from
        :py:func:`mpas_tools.viz.transect.horiz.make_triangle_tree()`

    subdivision_res : float, optional
        Resolution in m to use to subdivide the transect when looking for
        intersection candidates.  Should be small enough that curvature is
        small.

    Returns
    -------
    ds_out : xarray.Dataset
        A dataset that contains "nodes" where the transect intersects the
        edges of the triangles in ``ds_tris``.  The nodes also include the two
        end points of the transect, which typically lie within triangles. Each
        internal node (that is, not including the end points) is purposefully
        repeated twice, once for each triangle that node touches.  This allows
        for discontinuous fields between triangles (e.g. if one wishes to plot
        constant values on each MPAS cell).  The planar coordinates of these
        nodes are ``xNode`` and ``yNode``.  The distance along the transect of
        each intersection is ``dNode``. The index of the triangle and the first
        triangle node in ``ds_tris`` associated with each intersection node are
        given by ``horizTriangleIndices`` and ``horizTriangleNodeIndices``,
        respectively. The second node on the triangle for the edge associated
        with the intersection is given by
        ``numpy.mod(horizTriangleNodeIndices + 1, 3)``.

        The MPAS cell that a given node belongs to is given by
        ``horizCellIndices``. Each node also has an associated set of 6
        ``interpHorizCellIndices`` and ``interpHorizCellWeights`` that can be
        used to interpolate from MPAS cell centers to nodes first with
        area-weighted averaging to MPAS vertices and then linear interpolation
        along triangle edges.  Some of the weights may be zero, in which case
        the associated ``interpHorizCellIndices`` will be -1.

        Finally, ``xTransect`` and ``yTransect`` are included in the
        dataset, along with ``dTransect``, the distance along the transect of
        each original transect point.  In order to interpolate values (e.g.
        observations) from the original transect points to the intersection
        nodes, linear interpolation indices ``transectIndicesOnHorizNode`` and
        weights ``transectWeightsOnHorizNode`` are provided.  The values at
        nodes are found by::

          nodeValues = ((transectValues[transectIndicesOnHorizNode] *
                         transectWeightsOnHorizNode)
                        + (transectValues[transectIndicesOnHorizNode+1] *
                           (1.0 - transectWeightsOnHorizNode))
    """
    buffer = np.maximum(
        np.amax(ds_mesh.dvEdge.values), np.amax(ds_mesh.dcEdge.values)
    )

    n_nodes = ds_tris.sizes['nNodes']
    node_cell_weights = ds_tris.nodeCellWeights.values
    node_cell_indices = ds_tris.nodeCellIndices.values

    x = x_transect
    y = y_transect

    x_node = ds_tris.xNode.values.ravel()
    y_node = ds_tris.yNode.values.ravel()

    coordNode = np.zeros((len(x_node), 2))
    coordNode[:, 0] = x_node
    coordNode[:, 1] = y_node

    d_transect = np.zeros(x_transect.shape)

    d_node = None
    x_out = np.array([])
    y_out = np.array([])
    tris = None
    nodes = None
    interp_cells = None
    cell_weights = None

    n_horiz_weights = 6

    first = True

    d_start = 0.0
    for seg_index in range(len(x) - 1):
        sub_slice = slice(seg_index, seg_index + 2)
        x_sub, y_sub, _, _ = subdivide_planar(
            x[sub_slice].values, y[sub_slice].values, subdivision_res
        )

        start_point = Point(
            x_transect[seg_index].values, y_transect[seg_index].values
        )
        end_point = Point(
            x_transect[seg_index + 1].values, y_transect[seg_index + 1].values
        )

        segment = LineString([start_point, end_point])

        coords = np.zeros((len(x_sub), 3))
        coords[:, 0] = x_sub
        coords[:, 1] = y_sub
        radius = buffer + subdivision_res

        index_list = tree.query_ball_point(x=coords, r=radius)

        unique_indices = set()
        for indices in index_list:
            unique_indices.update(indices)

        start_indices = np.array(list(unique_indices))

        if len(start_indices) == 0:
            continue

        tris_cand = start_indices // n_nodes
        next_node_index = np.mod(start_indices + 1, n_nodes)
        end_indices = n_nodes * tris_cand + next_node_index

        intersecting_nodes = list()
        tris_inter_list = list()
        x_intersection_list = list()
        y_intersection_list = list()
        node_weights_list = list()
        node0_inter_list = list()
        node1_inter_list = list()
        distances_list = list()

        for index in range(len(start_indices)):
            start = start_indices[index]
            end = end_indices[index]

            node0 = Point(coordNode[start, 0], coordNode[start, 1])
            node1 = Point(coordNode[end, 0], coordNode[end, 1])

            edge = LineString([node0, node1])
            if segment.intersects(edge):
                point = segment.intersection(edge)
                intersecting_nodes.append((node0, node1, start, end, edge))

                if isinstance(point, LineString):
                    raise ValueError(
                        'A triangle edge exactly coincides with '
                        "a transect segment and I can't handle "
                        'that case.  Try moving the transect a '
                        'tiny bit.'
                    )
                elif not isinstance(point, Point):
                    raise ValueError(f'Unexpected intersection type {point}')

                x_intersection_list.append(point.x)
                y_intersection_list.append(point.y)

                start_to_intersection = LineString([start_point, point])

                weight = (
                    LineString([point, node1]).length
                    / LineString([node0, node1]).length
                )

                node_weights_list.append(weight)
                node0_inter_list.append(np.mod(start, n_nodes))
                node1_inter_list.append(np.mod(end, n_nodes))
                distances_list.append(start_to_intersection.length)
                tris_inter_list.append(tris_cand[index])

        distances = np.array(distances_list)
        x_intersection = np.array(x_intersection_list)
        y_intersection = np.array(y_intersection_list)
        node_weights = np.array(node_weights_list)
        node0_inter = np.array(node0_inter_list, dtype=int)
        node1_inter = np.array(node1_inter_list, dtype=int)
        tris_inter = np.array(tris_inter_list, dtype=int)

        d_node_local = d_start + distances

        d_start += segment.length

        weights = np.zeros((len(tris_inter), n_horiz_weights))
        cell_indices = np.zeros((len(tris_inter), n_horiz_weights), int)
        for index in range(3):
            weights[:, index] = (
                node_weights
                * node_cell_weights[tris_inter, node0_inter, index]
            )
            cell_indices[:, index] = node_cell_indices[
                tris_inter, node0_inter, index
            ]
            weights[:, index + 3] = (1.0 - node_weights) * node_cell_weights[
                tris_inter, node1_inter, index
            ]
            cell_indices[:, index + 3] = node_cell_indices[
                tris_inter, node1_inter, index
            ]

        if first:
            x_out = x_intersection
            y_out = y_intersection
            d_node = d_node_local

            tris = tris_inter
            nodes = node0_inter
            interp_cells = cell_indices
            cell_weights = weights
            first = False
        else:
            x_out = np.append(x_out, x_intersection)
            y_out = np.append(y_out, y_intersection)
            d_node = np.append(d_node, d_node_local)

            tris = np.concatenate((tris, tris_inter))
            nodes = np.concatenate((nodes, node0_inter))
            interp_cells = np.concatenate((interp_cells, cell_indices), axis=0)
            cell_weights = np.concatenate((cell_weights, weights), axis=0)

        d_transect[seg_index + 1] = d_start

    z_out = np.zeros(x_out.shape)

    epsilon = 1e-6 * subdivision_res
    (
        d_node,
        x_out,
        y_out,
        z_out,
        seg_tris,
        seg_nodes,
        interp_cells,
        cell_weights,
        valid_nodes,
    ) = _sort_intersections(
        d_node,
        tris,
        nodes,
        x_out,
        y_out,
        z_out,
        interp_cells,
        cell_weights,
        epsilon,
    )

    valid_segs = seg_tris >= 0
    cell_indices = -1 * np.ones(seg_tris.shape, dtype=int)
    cell_indices[valid_segs] = ds_tris.triCellIndices.values[
        seg_tris[valid_segs]
    ]

    ds_out = xr.Dataset()
    ds_out['xNode'] = (('nNodes',), x_out)
    ds_out['yNode'] = (('nNodes',), y_out)
    ds_out['dNode'] = (('nNodes',), d_node)

    ds_out['horizTriangleIndices'] = ('nSegments', seg_tris)
    ds_out['horizCellIndices'] = ('nSegments', cell_indices)
    ds_out['horizTriangleNodeIndices'] = (
        ('nSegments', 'nHorizBounds'),
        seg_nodes,
    )
    ds_out['interpHorizCellIndices'] = (
        ('nNodes', 'nHorizWeights'),
        interp_cells,
    )
    ds_out['interpHorizCellWeights'] = (
        ('nNodes', 'nHorizWeights'),
        cell_weights,
    )
    ds_out['validNodes'] = (('nNodes',), valid_nodes)

    transect_indices_on_horiz_node = np.zeros(d_node.shape, int)
    transect_weights_on_horiz_node = np.zeros(d_node.shape)
    for trans_index in range(len(d_transect) - 1):
        d0 = d_transect[trans_index]
        d1 = d_transect[trans_index + 1]
        mask = np.logical_and(d_node >= d0, d_node < d1)
        transect_indices_on_horiz_node[mask] = trans_index
        transect_weights_on_horiz_node[mask] = (d1 - d_node[mask]) / (d1 - d0)
    # last index will get missed by the mask and needs to be handled as a
    # special case
    transect_indices_on_horiz_node[-1] = len(d_transect) - 2
    transect_weights_on_horiz_node[-1] = 0.0

    ds_out['xTransect'] = x
    ds_out['yTransect'] = y
    ds_out['dTransect'] = (x_transect.dims, d_transect)
    ds_out['transectIndicesOnHorizNode'] = (
        ('nNodes',),
        transect_indices_on_horiz_node,
    )
    ds_out['transectWeightsOnHorizNode'] = (
        ('nNodes',),
        transect_weights_on_horiz_node,
    )

    return ds_out


def interp_mpas_horiz_to_transect_nodes(ds_transect, da):
    """
    Interpolate a 2D (``nCells``) MPAS DataArray to transect nodes, linearly
    interpolating fields between the closest neighboring cells

    Parameters
    ----------
    ds_transect : xr.Dataset
        A dataset that defines an MPAS transect, the results of calling
        ``find_spherical_transect_cells_and_weights()`` or
        ``find_planar_transect_cells_and_weights()``

    da : xr.DataArray
        An MPAS 2D field with dimensions `nCells`` (possibly among others)

    Returns
    -------
    da_nodes : xr.DataArray
        The data array interpolated to transect nodes with dimensions
        ``nNodes`` (in addition to whatever dimensions were in ``da`` besides
        ``nCells``)
    """
    interp_cell_indices = ds_transect.interpHorizCellIndices
    interp_cell_weights = ds_transect.interpHorizCellWeights
    da = da.isel(nCells=interp_cell_indices)
    da_nodes = (da * interp_cell_weights).sum(dim='nHorizWeights')

    da_nodes = da_nodes.where(ds_transect.validNodes)

    return da_nodes


def _sort_intersections(
    d_node,
    tris,
    nodes,
    x_out,
    y_out,
    z_out,
    interp_cells,
    cell_weights,
    epsilon,
):
    """sort nodes by distance and define segment between them"""

    sort_indices = np.argsort(d_node)
    d_sorted = d_node[sort_indices]

    # make a list of indices for each unique value of d
    d = d_sorted[0]
    unique_d_indices = [sort_indices[0]]
    unique_d_all_indices = [[sort_indices[0]]]
    for (
        index,
        next_d,
    ) in zip(sort_indices[1:], d_sorted[1:]):
        if next_d - d < epsilon:
            # this d value is effectively the same as the last, so we'll treat
            # it as the same
            unique_d_all_indices[-1].append(index)
        else:
            # this is a new d, so we'll add to a new list
            d = next_d
            unique_d_indices.append(index)
            unique_d_all_indices.append([index])

    # there is a segment between each unique d, though some are invalid (do
    # not correspond to a triangle)
    seg_tris_list = list()
    seg_nodes_list = list()

    index0 = unique_d_indices[0]
    indices0 = unique_d_all_indices[0]
    d0 = d_node[index0]

    indices = [index0]
    ds = [d0]
    for seg_index in range(len(unique_d_all_indices) - 1):
        indices1 = unique_d_all_indices[seg_index + 1]
        index1 = unique_d_indices[seg_index + 1]
        d1 = d_node[index1]

        # are there any triangles in common between this d value and the next?
        tris0 = tris[indices0]
        tris1 = tris[indices1]
        both = set(tris0).intersection(set(tris1))

        if len(both) > 0:
            tri = both.pop()
            seg_tris_list.append(tri)
            indices.append(index1)
            ds.append(d1)

            # the triangle nodes are the 2 corresponding to the same triangle
            # in the original list
            index0 = indices0[np.where(tris0 == tri)[0][0]]
            index1 = indices1[np.where(tris1 == tri)[0][0]]
            seg_nodes_list.append([nodes[index0], nodes[index1]])
        else:
            # this is an invalid segment so we need to insert and extra invalid
            # node to allow for proper masking
            seg_tris_list.extend([-1, -1])
            seg_nodes_list.extend([[-1, -1], [-1, -1]])
            indices.extend([index0, index1])
            ds.extend([0.5 * (d0 + d1), d1])

        index0 = index1
        indices0 = indices1
        d0 = d1

    indices = np.array(indices, dtype=int)
    d_node = np.array(ds, dtype=float)
    seg_tris = np.array(seg_tris_list, dtype=int)
    seg_nodes = np.array(seg_nodes_list, dtype=int)

    valid_nodes = np.ones(len(indices), dtype=bool)
    valid_nodes[1:-1] = np.logical_or(seg_tris[0:-1] >= 0, seg_tris[1:] > 0)

    x_out = x_out[indices]
    y_out = y_out[indices]
    z_out = z_out[indices]

    interp_cells = interp_cells[indices, :]
    cell_weights = cell_weights[indices, :]

    return (
        d_node,
        x_out,
        y_out,
        z_out,
        seg_tris,
        seg_nodes,
        interp_cells,
        cell_weights,
        valid_nodes,
    )


def _fix_periodic_tris(ds_tris, periodic_var, period):
    """
    make sure the given node coordinate on tris is within one period of the
    cell center
    """
    coord_node = ds_tris[periodic_var].values
    coord_cell = coord_node[:, 0]
    n_triangles = ds_tris.sizes['nTriangles']
    copy_pos = np.zeros(coord_cell.shape, dtype=bool)
    copy_neg = np.zeros(coord_cell.shape, dtype=bool)
    for i_node in [1, 2]:
        mask = coord_node[:, i_node] - coord_cell > 0.5 * period
        copy_pos = np.logical_or(copy_pos, mask)
        coord_node[:, i_node][mask] = coord_node[:, i_node][mask] - period
        mask = coord_node[:, i_node] - coord_cell < -0.5 * period
        copy_neg = np.logical_or(copy_neg, mask)
        coord_node[:, i_node][mask] = coord_node[:, i_node][mask] + period

    pos_indices = np.nonzero(copy_pos)[0]
    neg_indices = np.nonzero(copy_neg)[0]
    tri_indices = np.append(
        np.append(np.arange(0, n_triangles), pos_indices), neg_indices
    )

    ds_new = xr.Dataset(ds_tris)
    ds_new[periodic_var] = (('nTriangles', 'nNodes'), coord_node)
    ds_new = ds_new.isel(nTriangles=tri_indices)
    coord_node = ds_new[periodic_var].values

    pos_slice = slice(n_triangles, n_triangles + len(pos_indices))
    coord_node[pos_slice, :] = coord_node[pos_slice, :] + period
    neg_slice = slice(
        n_triangles + len(pos_indices),
        n_triangles + len(pos_indices) + len(neg_indices),
    )
    coord_node[neg_slice, :] = coord_node[neg_slice, :] - period
    ds_new[periodic_var] = (('nTriangles', 'nNodes'), coord_node)
    return ds_new
