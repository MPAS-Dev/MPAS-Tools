import argparse

import cmocean  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from geometric_features import FeatureCollection, read_feature_collection

from mpas_tools.ocean.viz.inset import add_inset
from mpas_tools.ocean.viz.transect.vert import (
    compute_transect,
    interp_mpas_to_transect_cells,
    interp_mpas_to_transect_nodes,
)
from mpas_tools.viz.colormaps import register_sci_viz_colormaps


def plot_transect(
    ds_transect,
    mpas_field=None,
    out_filename=None,
    ax=None,
    title=None,
    vmin=None,
    vmax=None,
    colorbar_label=None,
    cmap=None,
    figsize=(12, 6),
    dpi=200,
    method='flat',
    outline_color='black',
    ssh_color=None,
    seafloor_color=None,
    interface_color=None,
    cell_boundary_color=None,
    linewidth=1.0,
    color_start_and_end=False,
    start_color='red',
    end_color='green',
):
    """
    plot a transect showing the field on the MPAS-Ocean mesh and save to a file

    Parameters
    ----------
    ds_transect : xarray.Dataset
        A transect dataset from
        :py:func:`mpas_tools.ocean.viz.transect.vert.compute_transect()`

    mpas_field : xarray.DataArray
        The MPAS-Ocean 3D field to plot

    out_filename : str, optional
        The png file to write out to

    ax : matplotlib.axes.Axes
        Axes to plot to if making a multi-panel figure

    title : str
        The title of the plot

    vmin : float, optional
        The minimum values for the colorbar

    vmax : float, optional
        The maximum values for the colorbar

    colorbar_label : str, optional
        The colorbar label, or ``None`` if no colorbar is to be included.
        Use an empty string to display a colorbar without a label.

    cmap : str, optional
        The name of a colormap to use

    figsize : tuple, optional
        The size of the figure in inches

    dpi : int, optional
        The dots per inch of the image

    method : {'flat', 'bilinear'}, optional
        The type of interpolation to use in plots.  ``flat`` means constant
        values over each MPAS cell.  ``bilinear`` means smooth interpolation
        between horizontally between cell centers and vertical between the
        middle of layers.

    outline_color : str or None, optional
        The color to use to outline the transect or ``None`` for no outline

    ssh_color : str or None, optional
        The color to use to plot the SSH (sea surface height) or ``None`` if
        not plotting the SSH (except perhaps as part of the outline)

    seafloor_color : str or None, optional
        The color to use to plot the seafloor depth or ``None`` if not plotting
        the seafloor depth (except perhaps as part of the outline)

    interface_color : str or None, optional
        The color to use to plot interfaces between layers or ``None`` if
        not plotting the layer interfaces

    cell_boundary_color : str or None, optional
        The color to use to plot vertical boundaries between cells or ``None``
        if not plotting cell boundaries.  Typically, ``cell_boundary_color``
        will be used along with ``interface_color`` to outline cells both
        horizontally and vertically.

    linewidth : float, optional
        The width of outlines, interfaces and cell boundaries

    color_start_and_end : bool, optional
        Whether to color the left and right axes of the transect, which is
        useful if the transect is also being plotted in an inset or on top of
        a horizontal field

    start_color : str, optional
        The color of left axis marking the start of the transect if
        ``plot_start_end == True``

    end_color : str, optional
        The color of right axis marking the end of the transect if
        ``plot_start_end == True``
    """

    if ax is None and out_filename is None:
        raise ValueError('One of ax or out_filename must be supplied')

    create_fig = ax is None
    if create_fig:
        plt.figure(figsize=figsize)
        ax = plt.subplot(111)

    z = ds_transect.zTransectNode
    x = 1e-3 * ds_transect.dNode.broadcast_like(z)

    if mpas_field is not None:
        if method == 'flat':
            transect_field = interp_mpas_to_transect_cells(
                ds_transect, mpas_field
            )
            shading = 'flat'
        elif method == 'bilinear':
            transect_field = interp_mpas_to_transect_nodes(
                ds_transect, mpas_field
            )
            shading = 'gouraud'
        else:
            raise ValueError(f'Unsupported method: {method}')

        pc = ax.pcolormesh(
            x.values,
            z.values,
            transect_field.values,
            shading=shading,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            zorder=0,
        )
        ax.autoscale(tight=True)
        if colorbar_label is not None:
            plt.colorbar(
                pc, extend='both', shrink=0.7, ax=ax, label=colorbar_label
            )

    _plot_interfaces(
        ds_transect,
        ax,
        interface_color,
        cell_boundary_color,
        ssh_color,
        seafloor_color,
        color_start_and_end,
        start_color,
        end_color,
        linewidth,
    )

    _plot_outline(x, z, ds_transect.validCells, ax, outline_color, linewidth)

    ax.set_xlabel('transect distance (km)')
    ax.set_ylabel('z (m)')

    if create_fig:
        if title is not None:
            plt.title(title)
        plt.savefig(out_filename, dpi=dpi, bbox_inches='tight', pad_inches=0.2)
        plt.close()


def plot_feature_transects(
    fc,
    ds,
    ds_mesh=None,
    variable_list=None,
    cmap=None,
    flip=False,
    write_netcdf=False,
    method='flat',
    add_z=False,
):
    """
    Plot images of the given variables on the given transects.  One image
    named ``<transect_name>_<variable_name>.png`` will be produced in the
    current directory for each transect and variable

    Parameters
    ----------
    fc : geometric_features.FeatureCollection
        The transects to plot

    ds : xarray.Dataset
        The MPAS-Ocean dataset to plot

    ds_mesh : xarray.Dataset, optional
        The MPAS-Ocean mesh to use for plotting, the same as ``ds`` by default

    variable_list : list of str, optional
        The variables to plot

    cmap : str, optional
        The name of a colormap to use

    flip : book, optional
        Whether to flip the x axes of all transect plot

    write_netcdf : bool, optional
        Whether to write a NetCDF file for the transect in addition to the
        image

    method : {'flat', 'bilinear'}, optional
        The type of interpolation to use in plots.  ``flat`` means constant
        values over each MPAS cell.  ``bilinear`` means smooth interpolation
        between horizontally between cell centers and vertical between the
        middle of layers.

    add_z : bool, optional
        Whether to add zMid and zInterface to the mesh dataset
    """
    if 'Time' in ds.dims:
        ds = ds.isel(Time=0)

    if 'Time' in ds_mesh.dims:
        ds_mesh = ds_mesh.isel(Time=0)

    if add_z:
        _add_z(ds_mesh)

    print('\nBuilding transect geometry...')
    transects = _compute_feature_transects(fc, ds_mesh, flip)

    fc_transects = dict()
    for transect in fc.features:
        transect_name = transect['properties']['name']
        fc_transects[transect_name] = FeatureCollection(features=[transect])

    register_sci_viz_colormaps()

    if variable_list is None:
        variable_list = list()
        for var_name in ds.data_vars:
            var = ds[var_name]
            if 'nCells' in var.dims and (
                'nVertLevels' in var.dims or 'nVertLevelsP1' in var.dims
            ):
                variable_list.append(var_name)

    print('\nPlotting...')
    for var_name in variable_list:
        if var_name in ds:
            var = ds[var_name]
        elif var_name in ds_mesh:
            var = ds_mesh[var_name]
        else:
            raise ValueError(
                f'{var_name} not found in either the main or the '
                f'mesh dataset (if any)'
            )
        assert 'nCells' in var.dims and (
            'nVertLevels' in var.dims or 'nVertLevelsP1' in var.dims
        )
        for transect_name, ds_transect in transects.items():
            print(f'  {transect_name} {var_name}')
            _plot_feature_transect(
                ds_transect,
                var,
                var_name,
                transect_name,
                cmap,
                fc_transects[transect_name],
                write_netcdf,
                method,
            )


def plot_feature_transects_main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-g',
        '--geojson',
        dest='geojson_filename',
        required=True,
        help='A geojson file with transects to plot',
    )
    parser.add_argument(
        '-m',
        '--mesh',
        dest='mesh_filename',
        help='An MPAS-Ocean mesh file.  If not specified, the '
        'MPAS-Ocean data file must contain the mesh.',
    )
    parser.add_argument(
        '-f',
        '--file',
        dest='filename',
        required=True,
        help='An MPAS-Ocean data file',
    )
    parser.add_argument(
        '-v',
        '--variable_list',
        dest='variable_list',
        nargs='+',
        help='List of variables to plot.  All variables on '
        'cells in the data file is the default.',
    )
    parser.add_argument(
        '-c',
        '--colormap',
        dest='colormap',
        help='A colormap to use for the plots, default '
        'depends on the field name.',
    )
    parser.add_argument(
        '--flip',
        dest='flip',
        action='store_true',
        help='Flip the x axis for all transects',
    )
    parser.add_argument(
        '--write_netcdf',
        dest='write_netcdf',
        action='store_true',
        help='Whether to write a NetCDF file for the transect '
        'in addition to the image',
    )
    parser.add_argument(
        '--method',
        dest='method',
        default='flat',
        help='The type of interpolation to use in plots. '
        'Options are "flat" and "bilinear"',
    )
    parser.add_argument(
        '--add_z',
        dest='add_z',
        action='store_true',
        help='Whether to add zMid and zInterface to the mesh',
    )

    args = parser.parse_args()

    fc = read_feature_collection(args.geojson_filename)
    ds = xr.open_dataset(args.filename)
    if args.mesh_filename is not None:
        ds_mesh = xr.open_dataset(args.mesh_filename)
    else:
        ds_mesh = ds

    variable_list = args.variable_list

    if 'Time' in ds.dims:
        ds = ds.isel(Time=0)

    if 'Time' in ds_mesh.dims:
        ds_mesh = ds_mesh.isel(Time=0)

    plot_feature_transects(
        fc=fc,
        ds=ds,
        ds_mesh=ds_mesh,
        variable_list=variable_list,
        cmap=args.colormap,
        flip=args.flip,
        write_netcdf=args.write_netcdf,
        method=args.method,
        add_z=args.add_z,
    )


def _plot_interfaces(
    ds_transect,
    ax,
    interface_color,
    cell_boundary_color,
    ssh_color,
    seafloor_color,
    color_start_and_end,
    start_color,
    end_color,
    linewidth,
):
    if cell_boundary_color is not None:
        x_bnd = 1e-3 * ds_transect.dCellBoundary.values.T
        z_bnd = ds_transect.zCellBoundary.values.T
        ax.plot(
            x_bnd,
            z_bnd,
            color=cell_boundary_color,
            linewidth=linewidth,
            zorder=1,
        )

    if interface_color is not None:
        x_int = 1e-3 * ds_transect.dInterfaceSegment.values.T
        z_int = ds_transect.zInterfaceSegment.values.T
        ax.plot(
            x_int, z_int, color=interface_color, linewidth=linewidth, zorder=2
        )

    if ssh_color is not None:
        valid = ds_transect.validNodes.any(dim='nVertNodes')
        x_ssh = 1e-3 * ds_transect.dNode.values
        z_ssh = ds_transect.ssh.where(valid).values
        ax.plot(x_ssh, z_ssh, color=ssh_color, linewidth=linewidth, zorder=4)

    if seafloor_color is not None:
        valid = ds_transect.validNodes.any(dim='nVertNodes')
        x_floor = 1e-3 * ds_transect.dNode.values
        z_floor = ds_transect.zSeafloor.where(valid).values
        ax.plot(
            x_floor,
            z_floor,
            color=seafloor_color,
            linewidth=linewidth,
            zorder=5,
        )

    if color_start_and_end:
        ax.spines['left'].set_color(start_color)
        ax.spines['left'].set_linewidth(4 * linewidth)
        ax.spines['right'].set_color(end_color)
        ax.spines['right'].set_linewidth(4 * linewidth)


def _plot_outline(
    x, z, valid_cells, ax, outline_color, linewidth, epsilon=1e-6
):
    if outline_color is not None:
        # add a buffer of invalid values around the edge of the domain
        # and make copies of each node.  The validity of each copy of the node
        # corresponds to the validity of the adjacent cell
        valid = np.zeros(
            (2 * valid_cells.shape[0] + 2, 2 * valid_cells.shape[1] + 2),
            dtype=float,
        )
        z_buf = np.zeros(valid.shape, dtype=float)
        x_buf = np.zeros(valid.shape, dtype=float)

        valid_cells = valid_cells.astype(float)

        # each interior node get the value from its cell
        valid[1:-2:2, 1:-2:2] = valid_cells
        valid[2:-1:2, 1:-2:2] = valid_cells
        valid[1:-2:2, 2:-1:2] = valid_cells
        valid[2:-1:2, 2:-1:2] = valid_cells

        z_buf[:-1:2, :-1:2] = z.values
        z_buf[1::2, :-1:2] = z.values
        z_buf[:-1:2, 1::2] = z.values
        z_buf[1::2, 1::2] = z.values

        x_buf[:-1:2, :-1:2] = x.values
        x_buf[1::2, :-1:2] = x.values
        x_buf[:-1:2, 1::2] = x.values
        x_buf[1::2, 1::2] = x.values

        ax.contour(
            x_buf,
            z_buf,
            valid,
            levels=[1.0 - epsilon],
            colors=outline_color,
            linewidths=linewidth,
            zorder=3,
        )


def _compute_feature_transects(fc, ds_mesh, flip):
    """
    build a sequence of triangles showing the transect intersecting mpas cells
    """

    transects = dict()

    layer_thickness = ds_mesh.layerThickness
    bottom_depth = ds_mesh.bottomDepth
    max_level_cell = ds_mesh.maxLevelCell - 1
    if 'minLevelCell' in ds_mesh:
        min_level_cell = ds_mesh.minLevelCell - 1
    else:
        min_level_cell = xr.zeros_like(max_level_cell)

    spherical = ds_mesh.attrs['on_a_sphere'] == 'YES'

    for transect in fc.features:
        transect_name = transect['properties']['name']
        print(f'  {transect_name}')
        assert transect['geometry']['type'] == 'LineString'

        coordinates = transect['geometry']['coordinates']
        transect_lon, transect_lat = zip(*coordinates)
        transect_lon = np.array(transect_lon)
        transect_lat = np.array(transect_lat)
        if flip:
            transect_lon = transect_lon[::-1]
            transect_lat = transect_lat[::-1]
        transect_lon = xr.DataArray(data=transect_lon, dims=('nPoints',))
        transect_lat = xr.DataArray(data=transect_lat, dims=('nPoints',))

        ds_mpas_transect = compute_transect(
            x=transect_lon,
            y=transect_lat,
            ds_horiz_mesh=ds_mesh,
            layer_thickness=layer_thickness,
            bottom_depth=bottom_depth,
            min_level_cell=min_level_cell,
            max_level_cell=max_level_cell,
            spherical=spherical,
        )

        ds_mpas_transect.compute()
        transects[transect_name] = ds_mpas_transect

    return transects


def _plot_feature_transect(
    ds_transect,
    mpas_field,
    var_name,
    transect_name,
    cmap,
    fc,
    write_netcdf,
    method,
):
    """
    plot a transect showing the field on the MPAS-Ocean mesh and save to a file
    """
    transect_prefix = transect_name.replace(' ', '_')
    units = None
    if 'units' in mpas_field.attrs:
        units = mpas_field.attrs['units']

    colormaps = dict(
        temperature='cmo.thermal',
        salinity='cmo.haline',
        density='cmo.dense',
    )
    if cmap is None:
        for contains, map_name in colormaps.items():
            if contains in var_name.lower():
                cmap = map_name

    if units is not None:
        colorbar_label = f'{var_name} ({units})'
    else:
        colorbar_label = f'{var_name}'

    fig = plt.figure(figsize=(12, 6))
    ax = plt.gca()
    plot_transect(
        ds_transect=ds_transect,
        mpas_field=mpas_field,
        ax=ax,
        title=f'{var_name} through {transect_name}',
        colorbar_label=colorbar_label,
        cmap=cmap,
        color_start_and_end=True,
        method=method,
    )

    plt.tight_layout(pad=0.5, h_pad=0.5, rect=[0.0, 0.0, 1.0, 1.0])

    add_inset(fig, fc)
    plt.savefig(f'{transect_prefix}_{var_name}.png', dpi=200)
    plt.close()
    if write_netcdf:
        ds_transect.to_netcdf(f'{transect_prefix}_{var_name}.nc')


def _add_z(ds_mesh):
    """
    Add zMid and zInterface to ``ds_mesh``, useful for debugging
    """

    layer_thickness = ds_mesh.layerThickness
    bottom_depth = ds_mesh.bottomDepth
    max_level_cell = ds_mesh.maxLevelCell - 1
    if 'minLevelCell' in ds_mesh:
        min_level_cell = ds_mesh.minLevelCell - 1
    else:
        min_level_cell = xr.zeros_like(max_level_cell)

    n_vert_levels = layer_thickness.sizes['nVertLevels']

    vert_index = xr.DataArray.from_dict(
        {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)}
    )

    cell_mask = np.logical_and(
        vert_index >= min_level_cell, vert_index <= max_level_cell
    )
    layer_thickness = layer_thickness.where(cell_mask)

    thickness_sum = layer_thickness.sum(dim='nVertLevels')
    thickness_cum_sum = layer_thickness.cumsum(dim='nVertLevels')
    z_surface = -bottom_depth + thickness_sum

    z_layer_bot = z_surface - thickness_cum_sum

    z_interface_list = [z_surface]
    for z_index in range(n_vert_levels):
        z_interface_list.append(z_layer_bot.isel(nVertLevels=z_index))

    z_interface = xr.concat(z_interface_list, dim='nVertLevelsP1')

    vert_index = xr.DataArray.from_dict(
        {'dims': ('nVertLevelsP1',), 'data': np.arange(n_vert_levels + 1)}
    )
    interface_mask = np.logical_and(
        vert_index >= min_level_cell, vert_index <= max_level_cell + 1
    )

    z_interface = z_interface.where(interface_mask).transpose(
        'nCells', 'nVertLevelsP1'
    )

    z_mid = z_layer_bot + 0.5 * layer_thickness

    z_mid = z_mid.where(cell_mask).transpose('nCells', 'nVertLevels')

    ds_mesh.coords['zMid'] = z_mid
    ds_mesh.coords['zInterface'] = z_interface
