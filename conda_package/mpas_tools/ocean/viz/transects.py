#!/usr/bin/env python
import argparse

import cmocean
import numpy as np
import xarray as xr
from geometric_features import read_feature_collection, FeatureCollection
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpas_tools.viz import mesh_to_triangles
from mpas_tools.viz.transects import find_transect_cells_and_weights, \
    make_triangle_tree
from mpas_tools.ocean.transects import find_transect_levels_and_weights, \
    interp_mpas_to_transect_triangles, get_outline_segments

from mpas_tools.ocean.viz.inset import add_inset
from mpas_tools.viz.colormaps import register_sci_viz_colormaps


def plot_ocean_transects(fc, ds, ds_mesh=None, variable_list=None, cmap=None,
                         flip=False):
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
    """
    if 'Time' in ds.dims:
        ds = ds.isel(Time=0)

    if 'Time' in ds_mesh.dims:
        ds_mesh = ds_mesh.isel(Time=0)

    transects = _compute_transects(fc, ds_mesh, flip)

    print('\nBuilding transect geometry...')
    fc_transects = dict()
    for transect in fc.features:
        transect_name = transect['properties']['name']
        print(f'  {transect_name}')
        fc_transects[transect_name] = FeatureCollection(features=[transect])

    register_sci_viz_colormaps()

    if variable_list is None:
        variable_list = list()
        for var_name in ds.data_vars:
            var = ds[var_name]
            if 'nCells' in var.dims and 'nVertLevels' in var.dims:
                variable_list.append(var_name)

    print('\nPlotting...')
    for var_name in variable_list:
        var = ds[var_name]
        assert 'nCells' in var.dims and 'nVertLevels' in var.dims
        for transect_name, ds_transect in transects.items():
            print(f'  {transect_name} {var_name}')
            _plot_transect(ds_transect, var, var_name, transect_name, cmap,
                           fc_transects[transect_name])


def main():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--geojson', dest='geojson_filename',
                        required=True,
                        help='A geojson file with transects to plot')
    parser.add_argument('-m', '--mesh', dest='mesh_filename',
                        help='An MPAS-Ocean mesh file.  If not specified, the '
                             'MPAS-Ocean data file must contain the mesh.')
    parser.add_argument('-f', '--file', dest='filename', required=True,
                        help='An MPAS-Ocean data file')
    parser.add_argument("-v", "--variable_list", dest="variable_list",
                        nargs='+',
                        help="List of variables to plot.  All variables on "
                             "cells in the data file is the default.")
    parser.add_argument('-c', '--colormap', dest='colormap',
                        help='A colormap to use for the plots, default '
                             'depends on the field name.')
    parser.add_argument('--flip', dest='flip', action='store_true',
                        help='Flip the x axis for all transects')

    args = parser.parse_args()

    fc = read_feature_collection(args.geojson_filename)
    ds = xr.open_dataset(args.filename)
    if args.mesh_filename is not None:
        ds_mesh = xr.open_dataset(args.mesh_filename)
    else:
        ds_mesh = ds

    variable_list = args.variable_list

    plot_ocean_transects(fc=fc, ds=ds, ds_mesh=ds_mesh,
                         variable_list=variable_list, cmap=args.colormap,
                         flip=args.flip)


def _compute_transects(fc, ds_mesh, flip):
    """
    build a sequence of triangles showing the transect intersecting mpas cells
    """

    ds_tris = mesh_to_triangles(ds_mesh)

    triangle_tree = make_triangle_tree(ds_tris)

    transects = dict()

    for transect in fc.features:
        transect_name = transect['properties']['name']
        assert transect['geometry']['type'] == 'LineString'

        coordinates = transect['geometry']['coordinates']
        transect_lon, transect_lat = zip(*coordinates)
        transect_lon = np.array(transect_lon)
        transect_lat = np.array(transect_lat)
        if flip:
            transect_lon = transect_lon[::-1]
            transect_lat = transect_lat[::-1]
        transect_lon = xr.DataArray(data=transect_lon,
                                    dims=('nPoints',))
        transect_lat = xr.DataArray(data=transect_lat,
                                    dims=('nPoints',))

        ds_mpas_transect = find_transect_cells_and_weights(
            transect_lon, transect_lat, ds_tris, ds_mesh,
            triangle_tree, degrees=True)

        ds_mpas_transect = find_transect_levels_and_weights(
            ds_mpas_transect, ds_mesh.layerThickness,
            ds_mesh.bottomDepth, ds_mesh.maxLevelCell - 1)

        if 'landIceFraction' in ds_mesh:
            interp_cell_indices = ds_mpas_transect.interpHorizCellIndices
            interp_cell_weights = ds_mpas_transect.interpHorizCellWeights
            land_ice_fraction = ds_mesh.landIceFraction.isel(
                nCells=interp_cell_indices)
            land_ice_fraction = (land_ice_fraction * interp_cell_weights).sum(
                dim='nHorizWeights')
            ds_mpas_transect['landIceFraction'] = land_ice_fraction

        ds_mpas_transect['x'] = ds_mpas_transect.dNode.isel(
            nSegments=ds_mpas_transect.segmentIndices,
            nHorizBounds=ds_mpas_transect.nodeHorizBoundsIndices)

        ds_mpas_transect['z'] = ds_mpas_transect.zTransectNode

        ds_mpas_transect.compute()
        transects[transect_name] = ds_mpas_transect

    return transects


def _plot_transect(ds_transect, mpas_field, var_name, transect_name, cmap, fc):
    """
    plot a transect showing the field on the MPAS-Ocean mesh and save to a file
    """
    transect_prefix = transect_name.replace(' ', '_')
    transect_field = interp_mpas_to_transect_triangles(ds_transect,
                                                       mpas_field)
    units = None
    if 'units' in mpas_field.attrs:
        units = mpas_field.attrs['units']

    x_outline, z_outline = get_outline_segments(ds_transect)
    x_outline = 1e-3 * x_outline

    colormaps = dict(
        temperature='cmo.thermal',
        salinity='cmo.haline',
        density='cmo.dense',
    )
    if cmap is None:
        for contains, map_name in colormaps.items():
            if contains in var_name.lower():
                cmap = map_name

    tri_mask = np.logical_not(transect_field.notnull().values)
    # if any node of a triangle is masked, the triangle is masked
    # tri_mask = np.amax(tri_mask, axis=1)

    triangulation_args = _get_ds_triangulation_args(ds_transect)

    triangulation_args['mask'] = tri_mask

    tris = Triangulation(**triangulation_args)
    fig = plt.figure(figsize=(12, 6))
    ax = plt.gca()
    plt.tripcolor(tris, facecolors=transect_field.values, shading='flat',
                  cmap=cmap)
    plt.plot(x_outline, z_outline, 'k')
    if units is not None:
        colorbar_label = f'{var_name} ({units})'
    else:
        colorbar_label = f'{var_name}'
    plt.colorbar(label=colorbar_label)
    plt.title(f'{var_name} through {transect_name}')
    plt.xlabel('x (km)')
    plt.ylabel('z (m)')

    # make a red start axis and green end axis to correspond to the dots
    # in the inset
    ax.spines['left'].set_color('red')
    ax.spines['right'].set_color('green')
    ax.spines['left'].set_linewidth(4)
    ax.spines['right'].set_linewidth(4)

    plt.tight_layout(pad=0.5, h_pad=0.5, rect=[0.0, 0.0, 1.0, 1.0])

    add_inset(fig, fc)
    plt.savefig(f'{transect_prefix}_{var_name}.png', dpi=200)
    plt.close()


def _get_ds_triangulation_args(ds_transect):
    """
    get arguments for matplotlib Triangulation from triangulation dataset
    """

    n_transect_triangles = ds_transect.sizes['nTransectTriangles']
    d_node = ds_transect.dNode.isel(
        nSegments=ds_transect.segmentIndices,
        nHorizBounds=ds_transect.nodeHorizBoundsIndices)
    x = 1e-3 * d_node.values.ravel()

    z_transect_node = ds_transect.zTransectNode
    y = z_transect_node.values.ravel()

    tris = np.arange(3 * n_transect_triangles).reshape(
        (n_transect_triangles, 3))
    triangulation_args = dict(x=x, y=y, triangles=tris)

    return triangulation_args


if __name__ == '__main__':
    main()
