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
    z_transect=None,
):
    """
    build a sequence of quads showing the transect intersecting mpas cells.
    This can be used to plot transects of fields with dimensions ``nCells`` and
    either ``nVertLevels`` (levels) or ``nVertLevelsP1`` (interfaces).

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

    z_transect : xarray.DataArray, optional
        the z coordinate of the transect (1D or 2D).  If 2D, it must have the
        same along-transect dimension as ``x`` and ``y``

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

    ds_horiz_transect = _mask_land_segments(ds_horiz_transect, max_level_cell)

    ds_transect = find_transect_levels_and_weights(
        ds_horiz_transect=ds_horiz_transect,
        layer_thickness=layer_thickness,
        bottom_depth=bottom_depth,
        min_level_cell=min_level_cell,
        max_level_cell=max_level_cell,
        z_transect=z_transect,
    )

    ds_transect.compute()

    return ds_transect


def find_transect_levels_and_weights(
    ds_horiz_transect,
    layer_thickness,
    bottom_depth,
    min_level_cell,
    max_level_cell,
    z_transect=None,
):
    """
    Construct a vertical coordinate for a transect produced by
    :py:func:`mpas_tools.viz.transect.horiz.find_spherical_transect_cells_and_weights()`
    or :py:func:`mpas_tools.viz.transect.horiz.find_planar_transect_cells_and_weights()`.
    Also, compute interpolation weights such that observations at points on the
    original transect and with vertical coordinate ``transectZ`` can be
    bilinearly interpolated to the nodes of the transect.

    The transect is a slice through triangular prisms, one for each edge of
    each MPAS cell, with the cell center and the edge's two vertices as
    corners.  Layer interfaces are computed at cell centers from
    ``bottom_depth`` and ``layer_thickness``.  At a vertex, each interface
    lies at the kite-area-weighted average of that interface over the
    adjacent cells in which it is valid.  If this would put an interface
    below the one beneath it (as can happen where the same interface lies
    at very different depths in neighboring cells), the interfaces at that
    vertex are replaced by the closest non-increasing ones in the
    least-squares sense.  Interfaces vary linearly along each side of a
    triangle, so a transect node, where the transect crosses a side,
    takes a linear combination of the interfaces at the two ends of that
    side.  Levels that are valid in only one of the two segments adjacent
    to a node end in a vertical wall at that node.

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

    z_transect : xarray.DataArray, optional
        the z coordinate of the transect (1D or 2D).  If 2D, it must have the
        same along-transect dimension as ``ds_horiz_transect.dTransect``.

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

            - ``zTransectNode``: the elevation of each transect node, at
              layer interfaces (even ``nVertNodes`` indices) and the
              middles of layers (odd indices).  Where a level is not valid
              in either segment adjacent to a node, it has zero thickness
              at that node.
            - ``ssh``, ``zSeafloor``: the elevation of the top (sea surface
              or ice draft) and bottom (seafloor) of the valid layers in the
              segments adjacent to each node, the higher top and lower
              bottom where the two segments differ.  These are NaN at nodes
              between invalid segments.

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

    _add_vert_coord_and_interp_data(
        ds_transect_cells,
        layer_thickness,
        bottom_depth,
        min_level_cell,
        max_level_cell,
    )

    d_interface_seg, z_interface_seg = _get_interface_segments(
        ds_transect_cells.zTransectNode,
        ds_transect_cells.dNode,
        ds_transect_cells.validCells,
    )

    ds_transect_cells['dInterfaceSegment'] = d_interface_seg
    ds_transect_cells['zInterfaceSegment'] = z_interface_seg

    d_cell_boundary, z_cell_boundary = _get_cell_boundary_segments(
        ds_transect_cells.ssh,
        ds_transect_cells.zSeafloor,
        ds_transect_cells.dNode,
        ds_transect_cells.horizCellIndices,
    )

    ds_transect_cells['dCellBoundary'] = d_cell_boundary
    ds_transect_cells['zCellBoundary'] = z_cell_boundary

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

    if z_transect is not None:
        _add_vertical_interpolation_of_transect_points(
            ds_transect_cells, z_transect
        )

    return ds_transect_cells


def interp_mpas_to_transect_cells(ds_transect, da):
    """
    Interpolate an MPAS-Ocean DataArray with dimensions ``nCells`` by
    either ``nVertLevels`` (levels) or ``nVertLevelsP1`` (interfaces) to
    transect cells

    Parameters
    ----------
    ds_transect : xarray.Dataset
        A dataset that defines an MPAS-Ocean transect, the results of calling
        ``find_transect_levels_and_weights()``

    da : xarray.DataArray
        An MPAS-Ocean field with dimensions `nCells`` and either
        ``nVertLevels`` or ``nVertLevelsP1`` (possibly among others)

    Returns
    -------
    da_cells : xarray.DataArray
        The data array interpolated to transect cells with dimensions
        ``nSegments`` and ``nHalfLevels`` (in addition to whatever
        dimensions were in ``da`` besides ``nCells`` and either
        ``nVertLevels`` or ``nVertLevelsP1``)
    """

    cell_indices = ds_transect.cellIndices
    if 'nVertLevels' in da.dims:
        level_indices = ds_transect.levelIndices
        da_cells = da.isel(nCells=cell_indices, nVertLevels=level_indices)
    elif 'nVertLevelsP1' in da.dims:
        intreface_indices = ds_transect.interfaceIndices
        da_cells = da.isel(
            nCells=cell_indices, nVertLevelsP1=intreface_indices
        )
    else:
        raise ValueError(
            'da must have dimensions nCells and either '
            'nVertLevels or nVertLevelsP1'
        )

    da_cells = da_cells.where(ds_transect.validCells)

    return da_cells


def interp_mpas_to_transect_nodes(ds_transect, da):
    """
    Interpolate an MPAS-Ocean DataArray with dimensions ``nCells`` by
    either ``nVertLevels`` or ``nVertLevelsP1`` to transect nodes, linearly
    interpolating fields between the closest neighboring cells

    Parameters
    ----------
    ds_transect : xarray.Dataset
        A dataset that defines an MPAS-Ocean transect, the results of calling
        ``find_transect_levels_and_weights()``

    da : xarray.DataArray
        An MPAS-Ocean field with dimensions `nCells`` and either
        ``nVertLevels`` or ``nVertLevelsP1`` (possibly among others)

    Returns
    -------
    da_nodes : xarray.DataArray
        The data array interpolated to transect nodes with dimensions
        ``nHorizNodes`` and ``nVertNodes`` (in addition to whatever
        dimensions were in ``da`` besides ``nCells`` and either
        ``nVertLevels`` or ``nVertLevelsP1``)
    """

    interp_cell_indices = ds_transect.interpCellIndices

    if 'nVertLevels' in da.dims:
        interp_level_indices = ds_transect.interpLevelIndices
        interp_weights = ds_transect.interpLevelWeights

        da = da.isel(
            nCells=interp_cell_indices, nVertLevels=interp_level_indices
        )

    elif 'nVertLevelsP1' in da.dims:
        interp_interface_indices = ds_transect.interpInterfaceIndices
        interp_weights = ds_transect.interpInterfaceWeights

        da = da.isel(
            nCells=interp_cell_indices, nVertLevelsP1=interp_interface_indices
        )
    else:
        raise ValueError(
            'da must have dimensions nCells and either '
            'nVertLevels or nVertLevelsP1'
        )

    da_nodes = (da * interp_weights).sum(dim=('nHorizWeights', 'nVertWeights'))

    da_nodes = da_nodes.where(ds_transect.validNodes)

    return da_nodes


def interp_transect_grid_to_transect_nodes(ds_transect, da):
    """
    Interpolate a DataArray on the original transect grid to nodes on the
    MPAS-Ocean transect.

    Parameters
    ----------
    ds_transect : xarray.Dataset
        A dataset that defines an MPAS-Ocean transect, the
        results of calling ``find_transect_levels_and_weights()`` with the
        ``z_transect`` parameter.

    da : xarray.DataArray
        An field on the original transect (defined at vertical locations
        corresponding to ``z_transect``)

    Returns
    -------
    da_nodes : xarray.DataArray
        The data array interpolated to transect nodes with dimensions
        ``nHorizNodes`` and ``nVertNodes``
    """

    horiz_dim = ds_transect.dTransect.dims[0]
    z_transect = ds_transect.zTransect
    vert_dim = None
    for dim in z_transect.dims:
        if dim != horiz_dim:
            vert_dim = dim
            break

    assert vert_dim is not None

    horiz_indices = ds_transect.transectIndicesOnHorizNode
    horiz_weights = ds_transect.transectWeightsOnHorizNode

    vert_indices = ds_transect.transectInterpVertIndices
    vert_weights = ds_transect.transectInterpVertWeights

    kwargs00 = {horiz_dim: horiz_indices, vert_dim: vert_indices}
    kwargs01 = {horiz_dim: horiz_indices, vert_dim: vert_indices + 1}
    kwargs10 = {horiz_dim: horiz_indices + 1, vert_dim: vert_indices}
    kwargs11 = {horiz_dim: horiz_indices + 1, vert_dim: vert_indices + 1}

    da_nodes = (
        horiz_weights * vert_weights * da.isel(**kwargs00)
        + horiz_weights * (1.0 - vert_weights) * da.isel(**kwargs01)
        + (1.0 - horiz_weights) * vert_weights * da.isel(**kwargs10)
        + (1.0 - horiz_weights) * (1.0 - vert_weights) * da.isel(**kwargs11)
    )

    mask = np.logical_and(horiz_indices != -1, vert_indices != -1)

    da_nodes = da_nodes.where(mask)

    return da_nodes


def _mask_land_segments(ds_horiz_transect, max_level_cell):
    """
    Trim segments in land cells (``max_level_cell < 0``) from the ends of
    the transect and mark those in between as invalid, keeping a segment
    between each pair of adjacent nodes
    """
    cell_indices = ds_horiz_transect.horizCellIndices.values
    ocean = cell_indices >= 0
    ocean[ocean] = max_level_cell.values[cell_indices[ocean]] >= 0
    if not np.any(ocean):
        raise ValueError('The transect does not cross any ocean cells.')

    first = np.argmax(ocean)
    end = len(ocean) - np.argmax(ocean[::-1])
    ds_horiz_transect = ds_horiz_transect.isel(
        nSegments=slice(first, end), nNodes=slice(first, end + 1)
    )
    land = ~ocean[first:end]
    if np.any(land):
        ds_horiz_transect = ds_horiz_transect.copy()
        for var_name in [
            'horizCellIndices',
            'horizTriangleIndices',
            'horizTriangleNodeIndices',
        ]:
            var = ds_horiz_transect[var_name]
            values = var.values.copy()
            values[land, ...] = -1
            ds_horiz_transect[var_name] = (var.dims, values)
    return ds_horiz_transect


def _add_vert_coord_and_interp_data(
    ds_transect, layer_thickness, bottom_depth, min_level_cell, max_level_cell
):
    """
    Add the vertical coordinate of the transect and the indices and weights
    for interpolating MPAS-Ocean fields to transect nodes.

    Each MPAS cell is made up of triangular prisms, one per edge, with the
    cell center and the edge's two vertices as corners.  Layer interfaces
    are known at cell centers.  At a vertex, each interface lies at the
    kite-area-weighted average of that interface in the adjacent cells in
    which it is valid.  Along each side of a triangle, interfaces vary
    linearly, so the interfaces at a transect node, where the transect
    crosses a triangle side, are a linear combination of those at the two
    ends of that side.  Fields are interpolated to nodes in the same way.
    """
    n_horiz_nodes = ds_transect.sizes['nHorizNodes']
    n_segments = ds_transect.sizes['nSegments']
    n_vert_levels = layer_thickness.sizes['nVertLevels']

    # we assume below that there is a segment (whether valid or invalid)
    # connecting each pair of adjacent nodes
    assert n_horiz_nodes == n_segments + 1

    interp_horiz_cell_indices = ds_transect.interpHorizCellIndices
    cell_indices = ds_transect.horizCellIndices

    valid_cells = _get_cell_mask(
        cell_indices, min_level_cell, max_level_cell, n_vert_levels
    )
    valid_cells = valid_cells.transpose('nSegments', 'nVertLevels').values

    # the levels and interfaces that are valid in a segment adjacent to each
    # node
    valid_level_nodes = np.zeros((n_horiz_nodes, n_vert_levels), dtype=bool)
    valid_level_nodes[0:-1, :] = valid_cells
    valid_level_nodes[1:, :] = np.logical_or(
        valid_level_nodes[1:, :], valid_cells
    )
    valid_interface_nodes = np.zeros(
        (n_horiz_nodes, n_vert_levels + 1), dtype=bool
    )
    valid_interface_nodes[:, 0:-1] = valid_level_nodes
    valid_interface_nodes[:, 1:] = np.logical_or(
        valid_interface_nodes[:, 1:], valid_level_nodes
    )

    endpoints = _get_endpoint_cells_and_weights(
        interp_horiz_cell_indices.values,
        ds_transect.interpHorizCellWeights.values,
    )

    z_interface_cells, valid_level_cells, valid_interface_cells = (
        _get_cell_interfaces(
            endpoints['unique_cells'],
            layer_thickness,
            bottom_depth,
            min_level_cell,
            max_level_cell,
        )
    )

    # interfaces at the two ends of the triangle side each node lies on
    z_interface_ends, _ = _average_over_valid_cells(
        z_interface_cells, valid_interface_cells, endpoints
    )
    # a vertex average over different cells for different interfaces can
    # put an interface below the one beneath it
    z_interface_ends = z_interface_ends.reshape((-1, n_vert_levels + 1))
    z_interface_ends = _make_non_increasing(z_interface_ends)
    z_interface_ends = z_interface_ends.reshape(
        (n_horiz_nodes, -1, n_vert_levels + 1)
    )

    z_interface_nodes = _interp_ends_to_nodes(
        z_interface_ends, endpoints['side_weights']
    )

    z_interface_nodes, ssh, z_seafloor = _fill_interfaces(
        z_interface_nodes, valid_interface_nodes
    )

    z_half_interface = np.zeros((n_horiz_nodes, 2 * n_vert_levels + 1))
    z_half_interface[:, 0::2] = z_interface_nodes
    z_half_interface[:, 1::2] = 0.5 * (
        z_interface_nodes[:, 0:-1] + z_interface_nodes[:, 1:]
    )

    valid_transect_cells = np.zeros(
        (n_segments, 2 * n_vert_levels), dtype=bool
    )
    valid_transect_cells[:, 0::2] = valid_cells
    valid_transect_cells[:, 1::2] = valid_cells

    level_indices = np.zeros(2 * n_vert_levels, dtype=int)
    level_indices[0::2] = np.arange(n_vert_levels)
    level_indices[1::2] = np.arange(n_vert_levels)

    interface_indices = np.zeros(2 * n_vert_levels, dtype=int)
    interface_indices[0::2] = np.arange(n_vert_levels)
    interface_indices[1::2] = np.arange(1, n_vert_levels + 1)

    n_vert_nodes = 2 * n_vert_levels + 1
    n_vert_weights = 2

    interp_level_indices = -1 * np.ones(
        (n_vert_nodes, n_vert_weights), dtype=int
    )
    interp_level_indices[1:, 0] = level_indices
    interp_level_indices[0:-1, 1] = level_indices

    interp_interface_indices = -1 * np.ones(
        (n_vert_nodes, n_vert_weights), dtype=int
    )
    interp_interface_indices[1:, 0] = interface_indices
    interp_interface_indices[0:-1, 1] = interface_indices

    # horizontal weights for fields on levels and on interfaces
    _, horiz_level_weights = _average_over_valid_cells(
        None, valid_level_cells, endpoints
    )
    _, horiz_interface_weights = _average_over_valid_cells(
        None, valid_interface_cells, endpoints
    )

    # vertical weights for fields on levels: linear in z between the middles
    # of layers
    vert_level_weights = _get_vert_level_weights(
        z_interface_nodes, valid_level_nodes
    )
    interp_half_level_weights = _combine_horiz_and_vert_weights(
        horiz_level_weights, vert_level_weights, interp_level_indices
    )

    # vertical weights for fields on interfaces: the values at interfaces
    # and their mean at the middles of layers
    vert_interface_weights = np.zeros(
        (n_horiz_nodes, n_vert_nodes, n_vert_weights)
    )
    vert_interface_weights[:, 0::2, :] = (
        0.5 * valid_interface_nodes[:, :, np.newaxis]
    )
    vert_interface_weights[:, 0, :] = [0.0, 1.0] * valid_interface_nodes[
        :, 0:1
    ]
    vert_interface_weights[:, -1, :] = [1.0, 0.0] * valid_interface_nodes[
        :, -1:
    ]
    vert_interface_weights[:, 1::2, :] = (
        0.5 * valid_level_nodes[:, :, np.newaxis]
    )
    interp_half_interface_weights = _combine_horiz_and_vert_weights(
        horiz_interface_weights,
        vert_interface_weights,
        interp_interface_indices,
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

    ds_transect['zTransectNode'] = (
        ('nHorizNodes', 'nVertNodes'),
        z_half_interface,
    )
    ds_transect['ssh'] = (('nHorizNodes',), ssh)
    ds_transect['zSeafloor'] = (('nHorizNodes',), z_seafloor)

    ds_transect['cellIndices'] = cell_indices
    ds_transect['levelIndices'] = (('nHalfLevels',), level_indices)
    ds_transect['interfaceIndices'] = (('nHalfLevels',), interface_indices)
    ds_transect['validCells'] = (
        ('nSegments', 'nHalfLevels'),
        valid_transect_cells,
    )

    ds_transect['interpCellIndices'] = interp_horiz_cell_indices
    ds_transect['interpLevelIndices'] = (
        ('nVertNodes', 'nVertWeights'),
        interp_level_indices,
    )
    ds_transect['interpLevelWeights'] = (
        ('nHorizNodes', 'nVertNodes', 'nHorizWeights', 'nVertWeights'),
        interp_half_level_weights,
    )
    ds_transect['interpInterfaceIndices'] = (
        ('nVertNodes', 'nVertWeights'),
        interp_interface_indices,
    )
    ds_transect['interpInterfaceWeights'] = (
        ('nHorizNodes', 'nVertNodes', 'nHorizWeights', 'nVertWeights'),
        interp_half_interface_weights,
    )
    ds_transect['validNodes'] = (('nHorizNodes', 'nVertNodes'), valid_nodes)


def _get_endpoint_cells_and_weights(
    interp_horiz_cell_indices, interp_horiz_cell_weights
):
    """
    Split the horizontal interpolation cells and weights of each node into
    the two ends of the triangle side the node lies on.  Each end is either
    a cell center (one cell) or a vertex (the cells sharing that vertex).
    """
    n_nodes, n_horiz_weights = interp_horiz_cell_indices.shape
    if n_horiz_weights % 2 != 0:
        raise ValueError(
            f'Expected an even number of horizontal weights, one set for '
            f'each end of a triangle side, but found {n_horiz_weights}.'
        )
    n_cells_per_end = n_horiz_weights // 2
    shape = (n_nodes, 2, n_cells_per_end)
    cells = interp_horiz_cell_indices.reshape(shape)
    weights = interp_horiz_cell_weights.reshape(shape)

    # the weight of each end along the side of the triangle.  It includes
    # the weights of any missing cells (index -1) so that the ends' averages
    # can be taken over existing cells only
    side_weights = weights.sum(axis=2)

    cell_exists = cells >= 0
    unique_cells, unique_inverse = np.unique(
        np.where(cell_exists, cells, 0), return_inverse=True
    )
    unique_inverse = unique_inverse.reshape(shape)

    return dict(
        cells=cells,
        weights=np.where(cell_exists, weights, 0.0),
        side_weights=side_weights,
        cell_exists=cell_exists,
        unique_cells=unique_cells,
        unique_inverse=unique_inverse,
    )


def _get_cell_interfaces(
    cells, layer_thickness, bottom_depth, min_level_cell, max_level_cell
):
    """
    Get the elevations of layer interfaces at the centers of the given
    cells, building up from the seafloor as MPAS-Ocean does, along with
    which levels and interfaces are valid
    """
    n_vert_levels = layer_thickness.sizes['nVertLevels']
    cells = xr.DataArray(cells, dims=('nCells',))
    thickness = (
        layer_thickness.isel(nCells=cells)
        .transpose('nCells', 'nVertLevels')
        .values
    )
    bottom = bottom_depth.isel(nCells=cells).values
    min_level = min_level_cell.isel(nCells=cells).values
    max_level = max_level_cell.isel(nCells=cells).values

    level_indices = np.arange(n_vert_levels)
    valid_levels = np.logical_and(
        level_indices >= min_level[:, np.newaxis],
        level_indices <= max_level[:, np.newaxis],
    )
    interface_indices = np.arange(n_vert_levels + 1)
    valid_interfaces = np.logical_and(
        interface_indices >= min_level[:, np.newaxis],
        interface_indices <= max_level[:, np.newaxis] + 1,
    )
    valid_interfaces = np.logical_and(
        valid_interfaces, (max_level >= min_level)[:, np.newaxis]
    )

    thickness = np.where(valid_levels, thickness, 0.0)
    thickness_below = np.cumsum(thickness[:, ::-1], axis=1)[:, ::-1]
    z_interface = np.zeros((len(cells), n_vert_levels + 1))
    z_interface[:, 0:-1] = thickness_below
    z_interface = z_interface - bottom[:, np.newaxis]
    z_interface = np.where(valid_interfaces, z_interface, np.nan)

    return z_interface, valid_levels, valid_interfaces


def _average_over_valid_cells(field_cells, valid_cells, endpoints):
    """
    Average a field at each end of the triangle sides that nodes lie on,
    using only the cells in which the field is valid.  Returns the averages
    (if ``field_cells`` is not ``None``) and the weight of each cell in the
    linear interpolation to the node.
    """
    inverse = endpoints['unique_inverse']
    valid = np.logical_and(
        endpoints['cell_exists'][..., np.newaxis], valid_cells[inverse]
    )
    weights = endpoints['weights'][..., np.newaxis] * valid
    weight_sum = weights.sum(axis=2)
    has_valid = weight_sum > 0.0
    weights = weights / np.where(has_valid, weight_sum, 1.0)[:, :, np.newaxis]

    field_ends = None
    if field_cells is not None:
        field = np.where(valid, field_cells[inverse], 0.0)
        field_ends = (weights * field).sum(axis=2)
        field_ends = np.where(has_valid, field_ends, np.nan)

    # weights for interpolation to the node, flattened back to
    # nHorizWeights
    node_weights = endpoints['side_weights'][..., np.newaxis, np.newaxis] * (
        weights
    )
    n_nodes = node_weights.shape[0]
    node_weights = node_weights.reshape((n_nodes, -1, valid_cells.shape[1]))
    return field_ends, node_weights


def _interp_ends_to_nodes(field_ends, side_weights):
    """
    Interpolate linearly from the two ends of a triangle side to the node
    that lies on it.  The result is NaN where the field is not valid at an
    end with a nonzero weight.
    """
    use = side_weights > 0.0
    contrib = np.where(
        use[:, :, np.newaxis],
        side_weights[:, :, np.newaxis] * field_ends,
        0.0,
    )
    return contrib.sum(axis=1)


def _make_non_increasing(z):
    """
    Make each row of ``z`` non-increasing (skipping NaNs) with the closest
    such row in the least-squares sense
    """
    z_filled = _forward_fill(z)
    rising = np.any(z_filled[:, 1:] > z_filled[:, 0:-1], axis=1)
    for row in np.nonzero(rising)[0]:
        finite = np.isfinite(z[row, :])
        z[row, finite] = _pool_adjacent_violators(z[row, finite])
    return z


def _pool_adjacent_violators(values):
    """
    The non-increasing sequence closest to ``values`` in the least-squares
    sense, found by pooling adjacent values that increase into their mean
    """
    means = []
    counts = []
    for value in values:
        means.append(value)
        counts.append(1)
        while len(means) > 1 and means[-2] < means[-1]:
            count = counts[-2] + counts[-1]
            mean = (means[-2] * counts[-2] + means[-1] * counts[-1]) / count
            means.pop()
            counts.pop()
            means[-1] = mean
            counts[-1] = count
    return np.repeat(means, counts)


def _forward_fill(a):
    """
    Fill NaNs in each row with the closest preceding finite value
    """
    n_rows, n_cols = a.shape
    indices = np.where(np.isfinite(a), np.arange(n_cols)[np.newaxis, :], 0)
    np.maximum.accumulate(indices, axis=1, out=indices)
    return a[np.arange(n_rows)[:, np.newaxis], indices]


def _fill_interfaces(z_interface, valid_interfaces):
    """
    Give interfaces that are not valid at a node the elevation of the
    closest valid interface, so that invalid layers have zero thickness, and
    find the sea-surface height and seafloor elevation at each node: the top
    and bottom of the valid layers in the adjacent segments
    """
    n_nodes = z_interface.shape[0]
    any_valid = valid_interfaces.any(axis=1)
    # nodes between invalid segments keep whatever interfaces are defined
    mask = np.logical_or(valid_interfaces, ~any_valid[:, np.newaxis])
    z = np.where(mask, z_interface, np.nan)
    z = _forward_fill(z)
    z = _forward_fill(z[:, ::-1])[:, ::-1]
    z = np.where(np.isfinite(z), z, 0.0)

    top = np.argmax(valid_interfaces, axis=1)
    bottom = (
        valid_interfaces.shape[1]
        - 1
        - np.argmax(valid_interfaces[:, ::-1], axis=1)
    )
    nodes = np.arange(n_nodes)
    ssh = np.where(any_valid, z[nodes, top], np.nan)
    z_seafloor = np.where(any_valid, z[nodes, bottom], np.nan)
    return z, ssh, z_seafloor


def _get_vert_level_weights(z_interface, valid_levels):
    """
    Vertical weights for interpolating fields on levels to transect nodes:
    the level's value at the middle of a layer and linear interpolation
    between the middles of adjacent valid layers at an interface
    """
    n_nodes, n_vert_levels = valid_levels.shape
    n_vert_nodes = 2 * n_vert_levels + 1
    thickness = z_interface[:, 0:-1] - z_interface[:, 1:]

    weights = np.zeros((n_nodes, n_vert_nodes, 2))
    weights[:, 1::2, 0] = 0.5 * valid_levels
    weights[:, 1::2, 1] = 0.5 * valid_levels

    # the level above and below each interface
    above = np.zeros((n_nodes, n_vert_levels + 1), dtype=bool)
    above[:, 1:] = valid_levels
    below = np.zeros((n_nodes, n_vert_levels + 1), dtype=bool)
    below[:, 0:-1] = valid_levels
    thickness_above = np.zeros((n_nodes, n_vert_levels + 1))
    thickness_above[:, 1:] = thickness
    thickness_below = np.zeros((n_nodes, n_vert_levels + 1))
    thickness_below[:, 0:-1] = thickness

    both = np.logical_and(above, below)
    total = thickness_above + thickness_below
    # each middle is half a layer from the interface, so the weight of
    # each layer is proportional to the thickness of the other
    weight_above = np.where(
        total > 0.0, thickness_below / np.where(total > 0.0, total, 1.0), 0.5
    )
    weights[:, 0::2, 0] = np.where(both, weight_above, above)
    weights[:, 0::2, 1] = np.where(both, 1.0 - weight_above, below)
    return weights


def _combine_horiz_and_vert_weights(horiz_weights, vert_weights, indices):
    """
    Combine horizontal weights on each level or interface with vertical
    weights into weights with dimensions nHorizNodes, nVertNodes,
    nHorizWeights and nVertWeights
    """
    valid = indices >= 0
    horiz = horiz_weights[:, :, np.where(valid, indices, 0)]
    horiz = np.where(valid, horiz, 0.0)
    # nHorizNodes, nHorizWeights, nVertNodes, nVertWeights
    weights = horiz * vert_weights[:, np.newaxis, :, :]
    return weights.transpose((0, 2, 1, 3))


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


def _add_vertical_interpolation_of_transect_points(ds_transect, z_transect):
    d_transect = ds_transect.dTransect
    # make sure z_transect is 2D
    z_transect, _ = xr.broadcast(z_transect, d_transect)

    assert len(z_transect.dims) == 2

    horiz_dim = d_transect.dims[0]
    vert_dim = None
    for dim in z_transect.dims:
        if dim != horiz_dim:
            vert_dim = dim
            break

    assert vert_dim is not None

    nz_transect = z_transect.sizes[vert_dim]

    horiz_indices = ds_transect.transectIndicesOnHorizNode
    horiz_weights = ds_transect.transectWeightsOnHorizNode
    kwargs0 = {horiz_dim: horiz_indices}
    kwargs1 = {horiz_dim: horiz_indices + 1}
    z_transect_at_horiz_nodes = horiz_weights * z_transect.isel(**kwargs0) + (
        1.0 - horiz_weights
    ) * z_transect.isel(**kwargs1)

    z_transect_node = ds_transect.zTransectNode

    transect_interp_vert_indices = -1 * np.ones(
        z_transect_node.shape, dtype=np.int32
    )
    transect_interp_vert_weights = np.zeros(z_transect_node.shape)

    kwargs = {vert_dim: 0}
    z0 = z_transect_at_horiz_nodes.isel(**kwargs)
    for z_index in range(nz_transect - 1):
        kwargs = {vert_dim: z_index + 1}
        z1 = z_transect_at_horiz_nodes.isel(**kwargs)
        mask = np.logical_and(z_transect_node <= z0, z_transect_node > z1)
        mask = mask.values
        weights = (z1 - z_transect_node) / (z1 - z0)

        transect_interp_vert_indices[mask] = z_index
        transect_interp_vert_weights[mask] = weights.values[mask]
        z0 = z1

    ds_transect['transectInterpVertIndices'] = (
        z_transect_node.dims,
        transect_interp_vert_indices,
    )

    ds_transect['transectInterpVertWeights'] = (
        z_transect_node.dims,
        transect_interp_vert_weights,
    )

    ds_transect['zTransect'] = z_transect
