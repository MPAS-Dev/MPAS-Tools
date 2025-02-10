import numpy as np
import xarray as xr

from mpas_tools.viz.transect.horiz import (
    find_planar_transect_cells_and_weights,
    find_spherical_transect_cells_and_weights,
    make_triangle_tree,
    mesh_to_triangles,
)


def compute_transect(x, y, ds_horiz_mesh, layer_thickness, bottom_depth,
                     min_level_cell, max_level_cell, spherical=False,
                     z_transect=None):
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
            x, y, ds_tris, ds_horiz_mesh, triangle_tree, degrees=True)
    else:
        ds_horiz_transect = find_planar_transect_cells_and_weights(
            x, y, ds_tris, ds_horiz_mesh, triangle_tree)

    # mask horizontal transect to valid cells (max_level_cell >= 0)
    cell_indices = ds_horiz_transect.horizCellIndices
    seg_mask = max_level_cell.isel(nCells=cell_indices).values >= 0
    node_mask = np.zeros(ds_horiz_transect.sizes['nNodes'], dtype=bool)
    node_mask[0:-1] = seg_mask
    node_mask[1:] = np.logical_or(node_mask[1:], seg_mask)

    ds_horiz_transect = ds_horiz_transect.isel(nSegments=seg_mask,
                                               nNodes=node_mask)

    ds_transect = find_transect_levels_and_weights(
        ds_horiz_transect=ds_horiz_transect, layer_thickness=layer_thickness,
        bottom_depth=bottom_depth, min_level_cell=min_level_cell,
        max_level_cell=max_level_cell, z_transect=z_transect)

    ds_transect.compute()

    return ds_transect


def find_transect_levels_and_weights(ds_horiz_transect, layer_thickness,
                                     bottom_depth, min_level_cell,
                                     max_level_cell, z_transect=None):
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
        raise ValueError('Please select a single time level in layer '
                         'thickness.')

    ds_transect_cells = ds_horiz_transect.rename({'nNodes': 'nHorizNodes'})

    _add_vert_coord_and_interp_data(
        ds_transect_cells, layer_thickness, bottom_depth, min_level_cell,
        max_level_cell)

    d_interface_seg, z_interface_seg = _get_interface_segments(
        ds_transect_cells.zTransectNode, ds_transect_cells.dNode,
        ds_transect_cells.validCells)

    ds_transect_cells['dInterfaceSegment'] = d_interface_seg
    ds_transect_cells['zInterfaceSegment'] = z_interface_seg

    d_cell_boundary, z_cell_boundary = _get_cell_boundary_segments(
        ds_transect_cells.ssh, ds_transect_cells.zSeafloor,
        ds_transect_cells.dNode, ds_transect_cells.horizCellIndices)

    ds_transect_cells['dCellBoundary'] = d_cell_boundary
    ds_transect_cells['zCellBoundary'] = z_cell_boundary

    dims = ['nSegments', 'nHalfLevels', 'nHorizNodes', 'nVertNodes',
            'nInterfaceSegments', 'nCellBoundaries', 'nHorizBounds',
            'nVertBounds', 'nHorizWeights', 'nVertWeights']
    for dim in ds_transect_cells.dims:
        if dim not in dims:
            dims.insert(0, dim)
    ds_transect_cells = ds_transect_cells.transpose(*dims)

    if z_transect is not None:
        _add_vertical_interpolation_of_transect_points(
            ds_transect_cells, z_transect)

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
        da_cells = da.isel(nCells=cell_indices,
                           nVertLevelsP1=intreface_indices)
    else:
        raise ValueError('da must have dimensions nCells and either '
                         'nVertLevels or nVertLevelsP1')

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

        da = da.isel(nCells=interp_cell_indices,
                     nVertLevels=interp_level_indices)

    elif 'nVertLevelsP1' in da.dims:
        interp_interface_indices = ds_transect.interpInterfaceIndices
        interp_weights = ds_transect.interpInterfaceWeights

        da = da.isel(nCells=interp_cell_indices,
                     nVertLevelsP1=interp_interface_indices)
    else:
        raise ValueError('da must have dimensions nCells and either '
                         'nVertLevels or nVertLevelsP1')

    da_nodes = (da * interp_weights).sum(
        dim=('nHorizWeights', 'nVertWeights'))

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
·    """

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
    kwargs01 = {horiz_dim: horiz_indices, vert_dim: vert_indices+1}
    kwargs10 = {horiz_dim: horiz_indices + 1, vert_dim: vert_indices}
    kwargs11 = {horiz_dim: horiz_indices + 1, vert_dim: vert_indices+1}

    da_nodes = (
        horiz_weights * vert_weights * da.isel(**kwargs00) +
        horiz_weights * (1.0 - vert_weights) * da.isel(**kwargs01) +
        (1.0 - horiz_weights) * vert_weights * da.isel(**kwargs10) +
        (1.0 - horiz_weights) * (1.0 - vert_weights) * da.isel(**kwargs11))

    mask = np.logical_and(horiz_indices != -1, vert_indices != -1)

    da_nodes = da_nodes.where(mask)

    return da_nodes


def _get_horiz_at_interp(field, interp_horiz_cell_indices):
    field_interp = field.isel(nCells=interp_horiz_cell_indices)
    return field_interp


def _get_horiz_and_mask_at_interp(field, interp_horiz_cell_indices,
                                  cell_mask_interp):
    field_interp = _get_horiz_at_interp(field, interp_horiz_cell_indices)

    field_interp = field_interp.where(cell_mask_interp, 0.)
    return field_interp


def _get_interp_level_weights(valid_cells, n_horiz_nodes, n_vert_levels,
                              interp_horiz_cell_indices,
                              interp_horiz_cell_weights, cell_mask_interp):

    valid_nodes = np.zeros((n_horiz_nodes, n_vert_levels), dtype=bool)
    valid_nodes[0:-1, :] = valid_cells
    valid_nodes[1:, :] = np.logical_or(valid_nodes[1:, :], valid_cells)

    valid_nodes = xr.DataArray(dims=('nHorizNodes', 'nVertLevels'),
                               data=valid_nodes)

    interp_mask = np.logical_and(interp_horiz_cell_indices > 0,
                                 cell_mask_interp)

    interp_level_weights = interp_mask * interp_horiz_cell_weights
    weight_sum = interp_level_weights.sum(dim='nHorizWeights')

    valid_weights = valid_nodes.broadcast_like(interp_level_weights)
    interp_level_weights = \
        (interp_level_weights / weight_sum).where(valid_weights)

    return interp_level_weights


def _get_interp_interface_weights(valid_cells, n_horiz_nodes, n_vert_levels,
                                  interp_horiz_cell_indices,
                                  interp_horiz_cell_weights,
                                  cell_interface_interp):

    valid_nodes = np.zeros((n_horiz_nodes, n_vert_levels + 1), dtype=bool)
    valid_nodes[0:-1, 0:-1] = valid_cells
    valid_nodes[1:, 0:-1] = np.logical_or(valid_nodes[1:, 0:-1], valid_cells)
    valid_nodes[0:-1, 1:] = np.logical_or(valid_nodes[0:-1, 1:], valid_cells)
    valid_nodes[1:, 1:] = np.logical_or(valid_nodes[1:, 1:], valid_cells)

    valid_nodes = xr.DataArray(dims=('nHorizNodes', 'nVertLevelsP1'),
                               data=valid_nodes)

    interp_mask = np.logical_and(interp_horiz_cell_indices > 0,
                                 cell_interface_interp)

    interp_interface_weights = interp_mask * interp_horiz_cell_weights
    weight_sum = interp_interface_weights.sum(dim='nHorizWeights')

    valid_weights = valid_nodes.broadcast_like(interp_interface_weights)
    interp_interface_weights = \
        (interp_interface_weights / weight_sum).where(valid_weights)

    return interp_interface_weights


def _get_vert_coord_at_interp_cells(layer_thickness, bottom_depth,
                                    interp_horiz_cell_indices,
                                    cell_mask_interp):

    bottom_depth_interp = _get_horiz_at_interp(bottom_depth,
                                               interp_horiz_cell_indices)
    layer_thickness_interp = _get_horiz_and_mask_at_interp(
        layer_thickness, interp_horiz_cell_indices, cell_mask_interp)

    ssh_interp = (-bottom_depth_interp +
                  layer_thickness_interp.sum(dim='nVertLevels'))

    return ssh_interp, layer_thickness_interp


def _interp_horiz(field_at_interp, interp_weights):
    field_transect = (field_at_interp *
                      interp_weights).sum(dim='nHorizWeights')
    return field_transect


def _add_vert_coord_and_interp_data(ds_transect, layer_thickness, bottom_depth,
                                    min_level_cell, max_level_cell):
    n_horiz_nodes = ds_transect.sizes['nHorizNodes']
    n_segments = ds_transect.sizes['nSegments']
    n_vert_levels = layer_thickness.sizes['nVertLevels']

    # we assume below that there is a segment (whether valid or invalid)
    # connecting each pair of adjacent nodes
    assert n_horiz_nodes == n_segments + 1

    interp_horiz_cell_indices = ds_transect.interpHorizCellIndices
    interp_horiz_cell_weights = ds_transect.interpHorizCellWeights
    cell_indices = ds_transect.horizCellIndices

    cell_mask_interp = _get_cell_mask(interp_horiz_cell_indices,
                                      min_level_cell, max_level_cell,
                                      n_vert_levels)

    valid_cells = _get_cell_mask(cell_indices, min_level_cell, max_level_cell,
                                 n_vert_levels)

    valid_cells = valid_cells.transpose('nSegments', 'nVertLevels').values

    interp_level_weights = _get_interp_level_weights(
        valid_cells, n_horiz_nodes, n_vert_levels, interp_horiz_cell_indices,
        interp_horiz_cell_weights, cell_mask_interp)

    (ssh_at_interp, layer_thickness_at_interp) = \
        _get_vert_coord_at_interp_cells(layer_thickness, bottom_depth,
                                        interp_horiz_cell_indices,
                                        cell_mask_interp)

    ssh_transect = _interp_horiz(ssh_at_interp, interp_horiz_cell_weights)
    layer_thickness_transect = _interp_horiz(layer_thickness_at_interp,
                                             interp_level_weights)

    z_bot = ssh_transect - layer_thickness_transect.cumsum(dim='nVertLevels')
    z_mid = z_bot + 0.5 * layer_thickness_transect

    z_half_interfaces = [ssh_transect]
    for z_index in range(n_vert_levels):
        z_half_interfaces.extend([z_mid.isel(nVertLevels=z_index),
                                  z_bot.isel(nVertLevels=z_index)])

    z_half_interface = xr.concat(z_half_interfaces, dim='nVertNodes')
    z_half_interface = z_half_interface.transpose('nHorizNodes', 'nVertNodes')

    z_seafloor = ssh_transect - layer_thickness_transect.sum(
        dim='nVertLevels')

    valid_transect_cells = np.zeros((n_segments, 2 * n_vert_levels),
                                    dtype=bool)
    valid_transect_cells[:, 0::2] = valid_cells
    valid_transect_cells[:, 1::2] = valid_cells
    valid_transect_cells = xr.DataArray(dims=('nSegments', 'nHalfLevels'),
                                        data=valid_transect_cells)

    level_indices = np.zeros(2 * n_vert_levels, dtype=int)
    level_indices[0::2] = np.arange(n_vert_levels)
    level_indices[1::2] = np.arange(n_vert_levels)
    level_indices = xr.DataArray(dims=('nHalfLevels',), data=level_indices)

    interface_indices = np.zeros(2 * n_vert_levels, dtype=int)
    interface_indices[0::2] = np.arange(n_vert_levels)
    interface_indices[1::2] = np.arange(1, n_vert_levels + 1)
    interface_indices = xr.DataArray(dims=('nHalfLevels',),
                                     data=interface_indices)

    ds_transect['zTransectNode'] = z_half_interface

    ds_transect['ssh'] = ssh_transect
    ds_transect['zSeafloor'] = z_seafloor

    ds_transect['cellIndices'] = cell_indices
    ds_transect['levelIndices'] = level_indices
    ds_transect['interfaceIndices'] = interface_indices
    ds_transect['validCells'] = valid_transect_cells

    n_horiz_nodes = interp_horiz_cell_indices.sizes['nHorizNodes']
    n_vert_levels = layer_thickness.sizes['nVertLevels']
    n_vert_nodes = 2 * n_vert_levels + 1
    n_vert_weights = 2

    interp_level_indices = -1 * np.ones((n_vert_nodes, n_vert_weights),
                                        dtype=int)
    interp_level_indices[1:, 0] = level_indices.values
    interp_level_indices[0:-1, 1] = level_indices.values

    interp_level_indices = xr.DataArray(dims=('nVertNodes', 'nVertWeights'),
                                        data=interp_level_indices)

    interp_interface_indices = -1 * np.ones((n_vert_nodes, n_vert_weights),
                                            dtype=int)
    interp_interface_indices[1:, 0] = interface_indices.values
    interp_interface_indices[0:-1, 1] = interface_indices.values

    interp_interface_indices = xr.DataArray(
        dims=('nVertNodes', 'nVertWeights'),
        data=interp_interface_indices)

    half_level_thickness = 0.5 * layer_thickness.isel(
        nCells=interp_horiz_cell_indices, nVertLevels=interp_level_indices)
    half_level_thickness = half_level_thickness.where(
        interp_level_indices >= 0, other=0.)

    # vertical weights are proportional to the half-level thickness
    interp_half_level_weights = (
        half_level_thickness * interp_level_weights.isel(
            nVertLevels=interp_level_indices))

    weight_sum = interp_half_level_weights.sum(dim=('nHorizWeights',
                                               'nVertWeights'))
    out_mask = (weight_sum > 0.).broadcast_like(interp_half_level_weights)
    interp_half_level_weights = (interp_half_level_weights /
                                 weight_sum).where(out_mask)

    interface_mask_interp = _get_interface_mask(interp_horiz_cell_indices,
                                                min_level_cell, max_level_cell,
                                                n_vert_levels)

    interp_interface_weights = _get_interp_interface_weights(
        valid_cells, n_horiz_nodes, n_vert_levels, interp_horiz_cell_indices,
        interp_horiz_cell_weights, interface_mask_interp)

    # vertical weights are proportional to the half-level thickness
    interp_half_interface_weights = (
        half_level_thickness * interp_interface_weights.isel(
            nVertLevelsP1=interp_interface_indices))

    weight_sum = interp_half_interface_weights.sum(dim=('nHorizWeights',
                                                        'nVertWeights'))
    out_mask = (weight_sum > 0.).broadcast_like(interp_half_interface_weights)
    interp_half_interface_weights = (interp_half_interface_weights /
                                     weight_sum).where(out_mask)

    valid_nodes = np.zeros((n_horiz_nodes, n_vert_nodes), dtype=bool)
    valid_nodes[0:-1, 0:-1] = valid_transect_cells
    valid_nodes[1:, 0:-1] = np.logical_or(valid_nodes[1:, 0:-1],
                                          valid_transect_cells)
    valid_nodes[0:-1, 1:] = np.logical_or(valid_nodes[0:-1, 1:],
                                          valid_transect_cells)
    valid_nodes[1:, 1:] = np.logical_or(valid_nodes[1:, 1:],
                                        valid_transect_cells)

    valid_nodes = xr.DataArray(dims=('nHorizNodes', 'nVertNodes'),
                               data=valid_nodes)

    ds_transect['interpCellIndices'] = interp_horiz_cell_indices
    ds_transect['interpLevelIndices'] = interp_level_indices
    ds_transect['interpLevelWeights'] = interp_half_level_weights
    ds_transect['interpInterfaceIndices'] = interp_interface_indices
    ds_transect['interpInterfaceWeights'] = interp_half_interface_weights
    ds_transect['validNodes'] = valid_nodes


def _get_cell_mask(cell_indices, min_level_cell, max_level_cell,
                   n_vert_levels):
    level_indices = xr.DataArray(data=np.arange(n_vert_levels),
                                 dims='nVertLevels')
    min_level_cell = min_level_cell.isel(nCells=cell_indices)
    max_level_cell = max_level_cell.isel(nCells=cell_indices)

    cell_mask = np.logical_and(
        level_indices >= min_level_cell,
        level_indices <= max_level_cell)

    cell_mask = np.logical_and(cell_mask, cell_indices >= 0)

    return cell_mask


def _get_interface_mask(cell_indices, min_level_cell, max_level_cell,
                        n_vert_levels):
    interface_indices = xr.DataArray(data=np.arange(n_vert_levels + 1),
                                     dims='nVertLevelsP1')
    min_level_cell = min_level_cell.isel(nCells=cell_indices)
    max_level_cell = max_level_cell.isel(nCells=cell_indices)

    interface_mask = np.logical_and(
        interface_indices >= min_level_cell,
        interface_indices <= max_level_cell + 1)

    interface_mask = np.logical_and(interface_mask, cell_indices >= 0)

    return interface_mask


def _get_interface_segments(z_half_interface, d_node, valid_transect_cells):

    d = d_node.broadcast_like(z_half_interface)
    z_interface = z_half_interface.values[:, 0::2]
    d = d.values[:, 0::2]

    n_segments = valid_transect_cells.sizes['nSegments']
    n_half_levels = valid_transect_cells.sizes['nHalfLevels']
    n_vert_levels = n_half_levels // 2

    valid_segs = np.zeros((n_segments, n_vert_levels + 1), dtype=bool)
    valid_segs[:, 0:-1] = valid_transect_cells.values[:, 1::2]
    valid_segs[:, 1:] = np.logical_or(valid_segs[:, 1:],
                                      valid_transect_cells.values[:, 0::2])

    n_interface_segs = np.count_nonzero(valid_segs)

    d_seg = np.zeros((n_interface_segs, 2))
    z_seg = np.zeros((n_interface_segs, 2))
    d_seg[:, 0] = d[0:-1, :][valid_segs]
    d_seg[:, 1] = d[1:, :][valid_segs]
    z_seg[:, 0] = z_interface[0:-1, :][valid_segs]
    z_seg[:, 1] = z_interface[1:, :][valid_segs]

    d_seg = xr.DataArray(dims=('nInterfaceSegments', 'nHorizBounds'),
                         data=d_seg)

    z_seg = xr.DataArray(dims=('nInterfaceSegments', 'nHorizBounds'),
                         data=z_seg)

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


def _get_interp_indices_and_weights(layer_thickness, interp_cell_indices,
                                    interp_cell_weights, level_indices,
                                    valid_transect_cells):
    n_horiz_nodes = interp_cell_indices.sizes['nHorizNodes']
    n_vert_levels = layer_thickness.sizes['nVertLevels']
    n_vert_nodes = 2 * n_vert_levels + 1
    n_vert_weights = 2

    interp_level_indices = -1 * np.ones((n_vert_nodes, n_vert_weights),
                                        dtype=int)
    interp_level_indices[1:, 0] = level_indices.values
    interp_level_indices[0:-1, 1] = level_indices.values

    interp_level_indices = xr.DataArray(dims=('nVertNodes', 'nVertWeights'),
                                        data=interp_level_indices)

    half_level_thickness = 0.5 * layer_thickness.isel(
        nCells=interp_cell_indices, nVertLevels=interp_level_indices)
    half_level_thickness = half_level_thickness.where(
        interp_level_indices >= 0, other=0.)

    # vertical weights are proportional to the half-level thickness
    interp_cell_weights = half_level_thickness * interp_cell_weights.isel(
        nVertLevels=interp_level_indices)

    valid_nodes = np.zeros((n_horiz_nodes, n_vert_nodes), dtype=bool)
    valid_nodes[0:-1, 0:-1] = valid_transect_cells
    valid_nodes[1:, 0:-1] = np.logical_or(valid_nodes[1:, 0:-1],
                                          valid_transect_cells)
    valid_nodes[0:-1, 1:] = np.logical_or(valid_nodes[0:-1, 1:],
                                          valid_transect_cells)
    valid_nodes[1:, 1:] = np.logical_or(valid_nodes[1:, 1:],
                                        valid_transect_cells)

    valid_nodes = xr.DataArray(dims=('nHorizNodes', 'nVertNodes'),
                               data=valid_nodes)

    weight_sum = interp_cell_weights.sum(dim=('nHorizWeights', 'nVertWeights'))
    out_mask = (weight_sum > 0.).broadcast_like(interp_cell_weights)
    interp_cell_weights = (interp_cell_weights / weight_sum).where(out_mask)

    return interp_level_indices, interp_cell_weights, valid_nodes


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
    z_transect_at_horiz_nodes = (
        horiz_weights * z_transect.isel(**kwargs0) +
        (1.0 - horiz_weights) * z_transect.isel(**kwargs1))

    z_transect_node = ds_transect.zTransectNode

    transect_interp_vert_indices = -1*np.ones(z_transect_node.shape,
                                              dtype=np.int32)
    transect_interp_vert_weights = np.zeros(z_transect_node.shape)

    kwargs = {vert_dim: 0}
    z0 = z_transect_at_horiz_nodes.isel(**kwargs)
    for z_index in range(nz_transect-1):
        kwargs = {vert_dim: z_index + 1}
        z1 = z_transect_at_horiz_nodes.isel(**kwargs)
        mask = np.logical_and(z_transect_node <= z0, z_transect_node > z1)
        mask = mask.values
        weights = (z1 - z_transect_node)/(z1 - z0)

        transect_interp_vert_indices[mask] = z_index
        transect_interp_vert_weights[mask] = weights.values[mask]
        z0 = z1

    ds_transect['transectInterpVertIndices'] = (
        z_transect_node.dims, transect_interp_vert_indices)

    ds_transect['transectInterpVertWeights'] = (
        z_transect_node.dims, transect_interp_vert_weights)

    ds_transect['zTransect'] = z_transect
