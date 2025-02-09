import numpy as np
import xarray as xr

from mpas_tools.viz.transect.horiz import (
    find_planar_transect_cells_and_weights,
    find_spherical_transect_cells_and_weights,
    make_triangle_tree,
    mesh_to_triangles,
)


def compute_transect(
    x,
    y,
    ds_horiz_mesh,
    layer_thickness,
    bottom_depth,
    min_level_cell,
    max_level_cell,
    spherical=False,
):
    """
    build a sequence of quads showing the transect intersecting mpas cells.
    This can be used to plot transects of fields with dimensions ``nCells`` and
    ``nVertLevels``

    Parameters
    ----------
    x : xarray.DataArray
        The x or longitude coordinate of the transect

    y : xarray.DataArray
        The y or latitude coordinate of the transect

    ds_horiz_mesh : xarray.Dataset
        The horizontal MPAS mesh to use for plotting

    layer_thickness : xarray.DataArray
        The layer thickness at a particular instant in time.
        `layerThickness.isel(Time=tidx)` to select a particular time index
        `tidx` if the original data array contains `Time`.

    bottom_depth : xarray.DataArray
        the (positive down) depth of the seafloor on the MPAS mesh

    min_level_cell : xarray.DataArray
        the vertical zero-based index of the sea surface on the MPAS mesh

    max_level_cell : xarray.DataArray
        the vertical zero-based index of the bathymetry on the MPAS mesh

    spherical : bool, optional
        Whether the x and y coordinates are latitude and longitude in degrees

    Returns
    -------
    ds_transect : xarray.Dataset
        The transect dataset, see
        :py:func:`mpas_tools.ocean.viz.transect.vert.find_transect_levels_and_weights()`
        for details
    """  # noqa: E501

    ds_tris = mesh_to_triangles(ds_horiz_mesh)

    triangle_tree = make_triangle_tree(ds_tris)

    if spherical:
        ds_horiz_transect = find_spherical_transect_cells_and_weights(
            x, y, ds_tris, ds_horiz_mesh, triangle_tree, degrees=True
        )
    else:
        ds_horiz_transect = find_planar_transect_cells_and_weights(
            x, y, ds_tris, ds_horiz_mesh, triangle_tree
        )

    # mask horizontal transect to valid cells (max_level_cell >= 0)
    cell_indices = ds_horiz_transect.horizCellIndices
    seg_mask = max_level_cell.isel(nCells=cell_indices).values >= 0
    node_mask = np.zeros(ds_horiz_transect.sizes['nNodes'], dtype=bool)
    node_mask[0:-1] = seg_mask
    node_mask[1:] = np.logical_or(node_mask[1:], seg_mask)

    ds_horiz_transect = ds_horiz_transect.isel(
        nSegments=seg_mask, nNodes=node_mask
    )

    ds_transect = find_transect_levels_and_weights(
        ds_horiz_transect=ds_horiz_transect,
        layer_thickness=layer_thickness,
        bottom_depth=bottom_depth,
        min_level_cell=min_level_cell,
        max_level_cell=max_level_cell,
    )

    ds_transect.compute()

    return ds_transect


def find_transect_levels_and_weights(
    ds_horiz_transect,
    layer_thickness,
    bottom_depth,
    min_level_cell,
    max_level_cell,
):
    """
    Construct a vertical coordinate for a transect produced by
    :py:func:`mpas_tools.viz.transect.horiz.find_spherical_transect_cells_and_weights()`
    or :py:func:`mpas_tools.viz.transect.horiz.find_planar_transect_cells_and_weights()`.
    Also, compute interpolation weights such that observations at points on the
    original transect and with vertical coordinate ``transectZ`` can be
    bilinearly interpolated to the nodes of the transect.

    Parameters
    ----------
    ds_horiz_transect : xarray.Dataset
        A dataset that defines nodes of the transect

    layer_thickness : xarray.DataArray
        layer thicknesses on the MPAS mesh

    bottom_depth : xarray.DataArray
        the (positive down) depth of the seafloor on the MPAS mesh

    min_level_cell : xarray.DataArray
        the vertical zero-based index of the sea surface on the MPAS mesh

    max_level_cell : xarray.DataArray
        the vertical zero-based index of the bathymetry on the MPAS mesh

    Returns
    -------
    ds_transect : xarray.Dataset
        A dataset that contains nodes and cells that make up a 2D transect.

        There are ``nSegments`` horizontal and ``nHalfLevels`` vertical
        transect cells (quadrilaterals), bounded by ``nHorizNodes`` horizontal
        and ``nVertNodes`` vertical nodes (corners).

        In addition to the variables and coordinates in the input
        ``ds_transect``, the output dataset contains:

            - ``validCells``, ``validNodes``: which transect cells and nodes
              are valid (above the bathymetry and below the sea surface)

            - zTransectNode: the vertical height of each triangle node
            - ssh, zSeaFloor: the sea-surface height and sea-floor height at
              each node of each transect segment

            - ``cellIndices``: the MPAS-Ocean cell of a given transect segment
            - ``levelIndices``: the MPAS-Ocean vertical level of a given
              transect level

            - ``interpCellIndices``, ``interpLevelIndices``: the MPAS-Ocean
              cells and levels from which the value at a given transect cell is
              interpolated.  This can involve up to
              ``nHorizWeights * nVertWeights = 12`` different cells and levels.
            - interpCellWeights: the weight to multiply each field value by
              to perform interpolation to a transect cell.

            - ``dInterfaceSegment``, ``zInterfaceSegment`` - segments that can
              be used to plot the interfaces between MPAS-Ocean layers

            - ``dCellBoundary``, ``zCellBoundary`` - segments that can
              be used to plot the vertical boundaries between MPAS-Ocean cells

        Interpolation of a DataArray from MPAS cells and levels to transect
        cells can be performed with
        :py:func:`mpas_tools.ocean.viz.transect.vert.interp_mpas_to_transect_cells()`.
        Similarly, interpolation to transect nodes can be performed with
        :py:func:`mpas_tools.ocean.viz.transect.vert.interp_mpas_to_transect_nodes()`.
    """  # noqa: E501
    if 'Time' in layer_thickness.dims:
        raise ValueError(
            'Please select a single time level in layer thickness.'
        )

    ds_transect_cells = ds_horiz_transect.rename({'nNodes': 'nHorizNodes'})

    (
        z_half_interface,
        ssh,
        z_seafloor,
        interp_cell_indices,
        interp_cell_weights,
        valid_transect_cells,
        level_indices,
    ) = _get_vertical_coordinate(
        ds_transect_cells,
        layer_thickness,
        bottom_depth,
        min_level_cell,
        max_level_cell,
    )

    ds_transect_cells['zTransectNode'] = z_half_interface

    ds_transect_cells['ssh'] = ssh
    ds_transect_cells['zSeafloor'] = z_seafloor

    ds_transect_cells['cellIndices'] = ds_transect_cells.horizCellIndices
    ds_transect_cells['levelIndices'] = level_indices
    ds_transect_cells['validCells'] = valid_transect_cells

    d_interface_seg, z_interface_seg = _get_interface_segments(
        z_half_interface, ds_transect_cells.dNode, valid_transect_cells
    )

    ds_transect_cells['dInterfaceSegment'] = d_interface_seg
    ds_transect_cells['zInterfaceSegment'] = z_interface_seg

    d_cell_boundary, z_cell_boundary = _get_cell_boundary_segments(
        ssh,
        z_seafloor,
        ds_transect_cells.dNode,
        ds_transect_cells.horizCellIndices,
    )

    ds_transect_cells['dCellBoundary'] = d_cell_boundary
    ds_transect_cells['zCellBoundary'] = z_cell_boundary

    interp_level_indices, interp_cell_weights, valid_nodes = (
        _get_interp_indices_and_weights(
            layer_thickness,
            interp_cell_indices,
            interp_cell_weights,
            level_indices,
            valid_transect_cells,
        )
    )

    ds_transect_cells['interpCellIndices'] = interp_cell_indices
    ds_transect_cells['interpLevelIndices'] = interp_level_indices
    ds_transect_cells['interpCellWeights'] = interp_cell_weights
    ds_transect_cells['validNodes'] = valid_nodes

    dims = [
        'nSegments',
        'nHalfLevels',
        'nHorizNodes',
        'nVertNodes',
        'nInterfaceSegments',
        'nCellBoundaries',
        'nHorizBounds',
        'nVertBounds',
        'nHorizWeights',
        'nVertWeights',
    ]
    for dim in ds_transect_cells.dims:
        if dim not in dims:
            dims.insert(0, dim)
    ds_transect_cells = ds_transect_cells.transpose(*dims)

    return ds_transect_cells


def interp_mpas_to_transect_cells(ds_transect, da):
    """
    Interpolate an MPAS-Ocean DataArray with dimensions ``nCells`` by
    ``nVertLevels`` to transect cells

    Parameters
    ----------
    ds_transect : xarray.Dataset
        A dataset that defines an MPAS-Ocean transect, the results of calling
        ``find_transect_levels_and_weights()``

    da : xarray.DataArray
        An MPAS-Ocean field with dimensions `nCells`` and ``nVertLevels``
        (possibly among others)

    Returns
    -------
    da_cells : xarray.DataArray
        The data array interpolated to transect cells with dimensions
        ``nSegments`` and ``nHalfLevels`` (in addition to whatever
        dimensions were in ``da`` besides ``nCells`` and ``nVertLevels``)
    """

    cell_indices = ds_transect.cellIndices
    level_indices = ds_transect.levelIndices

    da_cells = da.isel(nCells=cell_indices, nVertLevels=level_indices)
    da_cells = da_cells.where(ds_transect.validCells)

    return da_cells


def interp_mpas_to_transect_nodes(ds_transect, da):
    """
    Interpolate an MPAS-Ocean DataArray with dimensions ``nCells`` by
    ``nVertLevels`` to transect nodes, linearly interpolating fields between
    the closest neighboring cells

    Parameters
    ----------
    ds_transect : xarray.Dataset
        A dataset that defines an MPAS-Ocean transect, the results of calling
        ``find_transect_levels_and_weights()``

    da : xarray.DataArray
        An MPAS-Ocean field with dimensions `nCells`` and ``nVertLevels``
        (possibly among others)

    Returns
    -------
    da_nodes : xarray.DataArray
        The data array interpolated to transect nodes with dimensions
        ``nHorizNodes`` and ``nVertNodes`` (in addition to whatever
        dimensions were in ``da`` besides ``nCells`` and ``nVertLevels``)
    """

    interp_cell_indices = ds_transect.interpCellIndices
    interp_level_indices = ds_transect.interpLevelIndices
    interp_cell_weights = ds_transect.interpCellWeights

    da = da.isel(nCells=interp_cell_indices, nVertLevels=interp_level_indices)

    da_nodes = (da * interp_cell_weights).sum(
        dim=('nHorizWeights', 'nVertWeights')
    )

    da_nodes = da_nodes.where(ds_transect.validNodes)

    return da_nodes


def _get_vertical_coordinate(
    ds_transect, layer_thickness, bottom_depth, min_level_cell, max_level_cell
):
    n_horiz_nodes = ds_transect.sizes['nHorizNodes']
    n_segments = ds_transect.sizes['nSegments']
    n_vert_levels = layer_thickness.sizes['nVertLevels']

    # we assume below that there is a segment (whether valid or invalid)
    # connecting each pair of adjacent nodes
    assert n_horiz_nodes == n_segments + 1

    interp_horiz_cell_indices = ds_transect.interpHorizCellIndices
    interp_horiz_cell_weights = ds_transect.interpHorizCellWeights

    bottom_depth_interp = bottom_depth.isel(nCells=interp_horiz_cell_indices)
    layer_thickness_interp = layer_thickness.isel(
        nCells=interp_horiz_cell_indices
    )

    cell_mask_interp = _get_cell_mask(
        interp_horiz_cell_indices,
        min_level_cell,
        max_level_cell,
        n_vert_levels,
    )
    layer_thickness_interp = layer_thickness_interp.where(
        cell_mask_interp, 0.0
    )

    ssh_interp = -bottom_depth_interp + layer_thickness_interp.sum(
        dim='nVertLevels'
    )

    interp_mask = np.logical_and(
        interp_horiz_cell_indices > 0, cell_mask_interp
    )

    interp_cell_weights = interp_mask * interp_horiz_cell_weights
    weight_sum = interp_cell_weights.sum(dim='nHorizWeights')

    cell_indices = ds_transect.horizCellIndices

    valid_cells = _get_cell_mask(
        cell_indices, min_level_cell, max_level_cell, n_vert_levels
    )

    valid_cells = valid_cells.transpose('nSegments', 'nVertLevels').values

    valid_nodes = np.zeros((n_horiz_nodes, n_vert_levels), dtype=bool)
    valid_nodes[0:-1, :] = valid_cells
    valid_nodes[1:, :] = np.logical_or(valid_nodes[1:, :], valid_cells)

    valid_nodes = xr.DataArray(
        dims=('nHorizNodes', 'nVertLevels'), data=valid_nodes
    )

    valid_weights = valid_nodes.broadcast_like(interp_cell_weights)
    interp_cell_weights = (interp_cell_weights / weight_sum).where(
        valid_weights
    )

    layer_thickness_transect = (
        layer_thickness_interp * interp_cell_weights
    ).sum(dim='nHorizWeights')

    interp_mask = max_level_cell.isel(nCells=interp_horiz_cell_indices) >= 0
    interp_horiz_cell_weights = interp_mask * interp_horiz_cell_weights
    weight_sum = interp_horiz_cell_weights.sum(dim='nHorizWeights')
    interp_horiz_cell_weights = (interp_horiz_cell_weights / weight_sum).where(
        interp_mask
    )

    ssh_transect = (ssh_interp * interp_horiz_cell_weights).sum(
        dim='nHorizWeights'
    )

    z_bot = ssh_transect - layer_thickness_transect.cumsum(dim='nVertLevels')
    z_mid = z_bot + 0.5 * layer_thickness_transect

    z_half_interfaces = [ssh_transect]
    for z_index in range(n_vert_levels):
        z_half_interfaces.extend(
            [z_mid.isel(nVertLevels=z_index), z_bot.isel(nVertLevels=z_index)]
        )

    z_half_interface = xr.concat(z_half_interfaces, dim='nVertNodes')
    z_half_interface = z_half_interface.transpose('nHorizNodes', 'nVertNodes')

    z_seafloor = ssh_transect - layer_thickness_transect.sum(dim='nVertLevels')

    valid_transect_cells = np.zeros(
        (n_segments, 2 * n_vert_levels), dtype=bool
    )
    valid_transect_cells[:, 0::2] = valid_cells
    valid_transect_cells[:, 1::2] = valid_cells
    valid_transect_cells = xr.DataArray(
        dims=('nSegments', 'nHalfLevels'), data=valid_transect_cells
    )

    level_indices = np.zeros(2 * n_vert_levels, dtype=int)
    level_indices[0::2] = np.arange(n_vert_levels)
    level_indices[1::2] = np.arange(n_vert_levels)
    level_indices = xr.DataArray(dims=('nHalfLevels',), data=level_indices)

    return (
        z_half_interface,
        ssh_transect,
        z_seafloor,
        interp_horiz_cell_indices,
        interp_cell_weights,
        valid_transect_cells,
        level_indices,
    )


def _get_cell_mask(
    cell_indices, min_level_cell, max_level_cell, n_vert_levels
):
    level_indices = xr.DataArray(
        data=np.arange(n_vert_levels), dims='nVertLevels'
    )
    min_level_cell = min_level_cell.isel(nCells=cell_indices)
    max_level_cell = max_level_cell.isel(nCells=cell_indices)

    cell_mask = np.logical_and(
        level_indices >= min_level_cell, level_indices <= max_level_cell
    )

    cell_mask = np.logical_and(cell_mask, cell_indices >= 0)

    return cell_mask


def _get_interface_segments(z_half_interface, d_node, valid_transect_cells):
    d = d_node.broadcast_like(z_half_interface)
    z_interface = z_half_interface.values[:, 0::2]
    d = d.values[:, 0::2]

    n_segments = valid_transect_cells.sizes['nSegments']
    n_half_levels = valid_transect_cells.sizes['nHalfLevels']
    n_vert_levels = n_half_levels // 2

    valid_segs = np.zeros((n_segments, n_vert_levels + 1), dtype=bool)
    valid_segs[:, 0:-1] = valid_transect_cells.values[:, 1::2]
    valid_segs[:, 1:] = np.logical_or(
        valid_segs[:, 1:], valid_transect_cells.values[:, 0::2]
    )

    n_interface_segs = np.count_nonzero(valid_segs)

    d_seg = np.zeros((n_interface_segs, 2))
    z_seg = np.zeros((n_interface_segs, 2))
    d_seg[:, 0] = d[0:-1, :][valid_segs]
    d_seg[:, 1] = d[1:, :][valid_segs]
    z_seg[:, 0] = z_interface[0:-1, :][valid_segs]
    z_seg[:, 1] = z_interface[1:, :][valid_segs]

    d_seg = xr.DataArray(
        dims=('nInterfaceSegments', 'nHorizBounds'), data=d_seg
    )

    z_seg = xr.DataArray(
        dims=('nInterfaceSegments', 'nHorizBounds'), data=z_seg
    )

    return d_seg, z_seg


def _get_cell_boundary_segments(ssh, z_seafloor, d_node, cell_indices):
    n_horiz_nodes = d_node.sizes['nHorizNodes']

    cell_boundary = np.ones(n_horiz_nodes, dtype=bool)
    cell_boundary[1:-1] = cell_indices.values[0:-1] != cell_indices.values[1:]

    n_cell_boundaries = np.count_nonzero(cell_boundary)

    d_seg = np.zeros((n_cell_boundaries, 2))
    z_seg = np.zeros((n_cell_boundaries, 2))
    d_seg[:, 0] = d_node.values[cell_boundary]
    d_seg[:, 1] = d_seg[:, 0]
    z_seg[:, 0] = ssh[cell_boundary]
    z_seg[:, 1] = z_seafloor[cell_boundary]

    d_seg = xr.DataArray(dims=('nCellBoundaries', 'nVertBounds'), data=d_seg)

    z_seg = xr.DataArray(dims=('nCellBoundaries', 'nVertBounds'), data=z_seg)

    return d_seg, z_seg


def _get_interp_indices_and_weights(
    layer_thickness,
    interp_cell_indices,
    interp_cell_weights,
    level_indices,
    valid_transect_cells,
):
    n_horiz_nodes = interp_cell_indices.sizes['nHorizNodes']
    n_vert_levels = layer_thickness.sizes['nVertLevels']
    n_vert_nodes = 2 * n_vert_levels + 1
    n_vert_weights = 2

    interp_level_indices = -1 * np.ones(
        (n_vert_nodes, n_vert_weights), dtype=int
    )
    interp_level_indices[1:, 0] = level_indices.values
    interp_level_indices[0:-1, 1] = level_indices.values

    interp_level_indices = xr.DataArray(
        dims=('nVertNodes', 'nVertWeights'), data=interp_level_indices
    )

    half_level_thickness = 0.5 * layer_thickness.isel(
        nCells=interp_cell_indices, nVertLevels=interp_level_indices
    )
    half_level_thickness = half_level_thickness.where(
        interp_level_indices >= 0, other=0.0
    )

    # vertical weights are proportional to the half-level thickness
    interp_cell_weights = half_level_thickness * interp_cell_weights.isel(
        nVertLevels=interp_level_indices
    )

    valid_nodes = np.zeros((n_horiz_nodes, n_vert_nodes), dtype=bool)
    valid_nodes[0:-1, 0:-1] = valid_transect_cells
    valid_nodes[1:, 0:-1] = np.logical_or(
        valid_nodes[1:, 0:-1], valid_transect_cells
    )
    valid_nodes[0:-1, 1:] = np.logical_or(
        valid_nodes[0:-1, 1:], valid_transect_cells
    )
    valid_nodes[1:, 1:] = np.logical_or(
        valid_nodes[1:, 1:], valid_transect_cells
    )

    valid_nodes = xr.DataArray(
        dims=('nHorizNodes', 'nVertNodes'), data=valid_nodes
    )

    weight_sum = interp_cell_weights.sum(dim=('nHorizWeights', 'nVertWeights'))
    out_mask = (weight_sum > 0.0).broadcast_like(interp_cell_weights)
    interp_cell_weights = (interp_cell_weights / weight_sum).where(out_mask)

    return interp_level_indices, interp_cell_weights, valid_nodes
