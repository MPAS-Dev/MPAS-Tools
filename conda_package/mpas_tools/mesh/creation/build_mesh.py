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
from mpas_tools.viz.paraview_extractor import extract_vtk

from mpas_tools.mesh.creation.jigsaw_driver import jigsaw_driver
from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.viz.colormaps import register_sci_viz_colormaps


def build_spherical_mesh(cellWidth, lon, lat, earth_radius,
                         out_filename='base_mesh.nc', plot_cellWidth=True,
                         vtk_dir='base_mesh_vtk'):
    """
    Build an MPAS mesh using JIGSAW with the given cell sizes as a function of
    latitude and longitude (on a sphere) or x and y (on a plane).

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

    vtk_dir : str, optional
        The name of the directory where mesh data will be extracted for viewing
        in ParaVeiw.
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

    _shared_steps(cw_filename, out_filename, vtk_dir, on_sphere=True)


def build_planar_mesh(cellWidth, x, y, geom_points, geom_edges,
                      out_filename='base_mesh.nc', vtk_dir='base_mesh_vtk'):
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

    vtk_dir : str, optional
        The name of the directory where mesh data will be extracted for viewing
        in ParaVeiw.

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

    _shared_steps(cw_filename, out_filename, vtk_dir, on_sphere=False)


def _shared_steps(cw_filename, out_filename, vtk_dir, on_sphere):
    print('Step 2. Convert triangles from jigsaw format to netcdf')
    jigsaw_to_netcdf(msh_filename='mesh-MESH.msh',
                     output_name='mesh_triangles.nc', on_sphere=on_sphere)

    print('Step 3. Convert from triangles to MPAS mesh')
    write_netcdf(convert(xarray.open_dataset('mesh_triangles.nc')),
                 out_filename)

    print('Step 4. Create vtk file for visualization')
    extract_vtk(ignore_time=True, lonlat=True, dimension_list=['maxEdges='],
                variable_list=['allOnCells'], filename_pattern=out_filename,
                out_dir=vtk_dir)

    print("***********************************************")
    print("**    The global mesh file is {}   **".format(out_filename))
    print("***********************************************")
