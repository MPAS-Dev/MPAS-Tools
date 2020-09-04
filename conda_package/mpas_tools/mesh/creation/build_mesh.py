from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy

from mpas_tools.mesh.conversion import convert
from mpas_tools.io import write_netcdf

from mpas_tools.mesh.creation.jigsaw_driver import jigsaw_driver
from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.viz.colormaps import register_sci_viz_colormaps


def build_spherical_mesh(cellWidth, lon, lat, earth_radius,
                         out_filename='base_mesh.nc', plot_cellWidth=True):
    """
    Build an MPAS mesh using JIGSAW with the given cell sizes as a function of
    latitude and longitude.

    The result is a mesh file stored in ``out_filename`` as well as several
    intermediate files: ``mesh.log``, ``mesh-HFUN.msh``, ``mesh.jig``,
    ``mesh-MESH.msh``, ``mesh.msh``, and ``mesh_triangles.nc``.

    Parameters
    ----------
    cellWidth : ndarray
        m x n array of cell width in km

    lon : ndarray
        longitude in degrees (length n and between -180 and 180)

    lat : ndarray
        longitude in degrees (length m and between -90 and 90)

    earth_radius : float
        Earth radius in meters

    out_filename : str, optional
        The file name of the resulting MPAS mesh

    plot_cellWidth : bool, optional
        Whether to produce a plot of ``cellWidth``. If so, it will be written to
        ``cellWidthGlobal.png``.
    """

    da = xarray.DataArray(cellWidth,
                          dims=['lat', 'lon'],
                          coords={'lat': lat, 'lon': lon},
                          name='cellWidth')
    cw_filename = 'cellWidthVsLatLon.nc'
    da.to_netcdf(cw_filename)
    if plot_cellWidth:
        register_sci_viz_colormaps()
        fig = plt.figure(figsize=[16.0, 8.0])
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_global()
        im = ax.imshow(cellWidth, origin='lower',
                       transform=ccrs.PlateCarree(),
                       extent=[-180, 180, -90, 90], cmap='3Wbgy5',
                       zorder=0)
        ax.add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)
        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=1,
            color='gray',
            alpha=0.5,
            linestyle='-', zorder=2)
        gl.top_labels = False
        gl.right_labels = False
        plt.title(
            'Grid cell size, km, min: {:.1f} max: {:.1f}'.format(
            cellWidth.min(),cellWidth.max()))
        plt.colorbar(im, shrink=.60)
        fig.canvas.draw()
        plt.tight_layout()
        plt.savefig('cellWidthGlobal.png', bbox_inches='tight')
        plt.close()

    print('Step 1. Generate mesh with JIGSAW')
    jigsaw_driver(cellWidth, lon, lat, on_sphere=True,
                  earth_radius=earth_radius)

    print('Step 2. Convert triangles from jigsaw format to netcdf')
    jigsaw_to_netcdf(msh_filename='mesh-MESH.msh',
                     output_name='mesh_triangles.nc', on_sphere=True,
                     sphere_radius=earth_radius)

    print('Step 3. Convert from triangles to MPAS mesh')
    write_netcdf(convert(xarray.open_dataset('mesh_triangles.nc')),
                 out_filename)


def build_planar_mesh(cellWidth, x, y, geom_points, geom_edges,
                      out_filename='base_mesh.nc'):
    """
    Build a planar MPAS mesh

    Parameters
    ----------
    cellWidth : ndarray
        m x n array of cell width in km

    x, y : ndarray
        arrays defining planar coordinates in meters

    geom_points : ndarray
        list of point coordinates for bounding polygon for the planar mesh

    geom_edges : ndarray
        list of edges between points in ``geom_points`` that define the
        bounding polygon

    out_filename : str, optional
        The file name of the resulting MPAS mesh
    """
    da = xarray.DataArray(cellWidth,
                          dims=['y', 'x'],
                          coords={'y': y, 'x': x},
                          name='cellWidth')
    cw_filename = 'cellWidthVsXY.nc'
    da.to_netcdf(cw_filename)

    print('Step 1. Generate mesh with JIGSAW')
    jigsaw_driver(
        cellWidth,
        x,
        y,
        on_sphere=False,
        geom_points=geom_points,
        geom_edges=geom_edges)

    print('Step 2. Convert triangles from jigsaw format to netcdf')
    jigsaw_to_netcdf(msh_filename='mesh-MESH.msh',
                     output_name='mesh_triangles.nc', on_sphere=False)

    print('Step 3. Convert from triangles to MPAS mesh')
    write_netcdf(convert(xarray.open_dataset('mesh_triangles.nc')),
                 out_filename)

