from __future__ import absolute_import, division, print_function, \
    unicode_literals

from mpas_tools.mesh.creation import build_spherical_mesh as create_spherical_mesh
from mpas_tools.mesh.creation import build_planar_mesh as create_planar_mesh
from mpas_tools.ocean.inject_bathymetry import inject_bathymetry
from mpas_tools.ocean.inject_meshDensity import inject_spherical_meshDensity, \
    inject_planar_meshDensity
from mpas_tools.ocean.inject_preserve_floodplain import inject_preserve_floodplain
from mpas_tools.viz.paraview_extractor import extract_vtk


def build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc',
                         plot_cellWidth=True, vtk_dir='base_mesh_vtk',
                         preserve_floodplain=False, floodplain_elevation=20.0,
                         do_inject_bathymetry=False):
    """
    Build an MPAS mesh using JIGSAW with the given cell sizes as a function of
    latitude and longitude

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

    out_filename : str, optional
        The file name of the resulting MPAS mesh

    plot_cellWidth : bool, optional
        Whether to produce a plot of ``cellWidth``. If so, it will be written to
        ``cellWidthGlobal.png``.

    vtk_dir : str, optional
        The name of the directory where mesh data will be extracted for viewing
        in ParaVeiw.

    preserve_floodplain : bool, optional
        Whether a flood plain (bathymetry above z = 0) should  be preserved in
        the mesh

    floodplain_elevation : float, optional
        The elevation in meters to which the flood plain is preserved

    do_inject_bathymetry : bool, optional
        Whether one of the default bathymetry datasets, ``earth_relief_15s.nc``
        or ``topo.msh``, should be added to the MPAS mesh
    """

    # from https://github.com/E3SM-Project/E3SM/blob/master/cime/src/share/util/shr_const_mod.F90#L20
    earth_radius = 6.37122e6

    create_spherical_mesh(cellWidth, lon, lat, earth_radius,  out_filename,
                          plot_cellWidth)

    print('Step 4. Inject meshDensity into the mesh file')
    inject_spherical_meshDensity(cellWidth, lon, lat,
                                 mesh_filename=out_filename)

    _shared_steps(out_filename, vtk_dir, preserve_floodplain,
                  floodplain_elevation, do_inject_bathymetry)


def build_planar_mesh(cellWidth, x, y, geom_points, geom_edges,
                      out_filename='base_mesh.nc', vtk_dir='base_mesh_vtk',
                      preserve_floodplain=False, floodplain_elevation=20.0,
                      do_inject_bathymetry=False):
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

    preserve_floodplain : bool, optional
        Whether a flood plain (bathymetry above z = 0) should  be preserved in
        the mesh

    floodplain_elevation : float, optional
        The elevation in meters to which the flood plain is preserved

    do_inject_bathymetry : bool, optional
        Whether one of the default bathymetry datasets, ``earth_relief_15s.nc``
        or ``topo.msh``, should be added to the MPAS mesh
    """
    create_planar_mesh(cellWidth, x, y, geom_points, geom_edges, out_filename)

    print('Step 4. Inject meshDensity into the mesh file')
    inject_planar_meshDensity(cellWidth, x, y, mesh_filename=out_filename)

    _shared_steps(out_filename, vtk_dir, preserve_floodplain,
                  floodplain_elevation, do_inject_bathymetry)


def _shared_steps(out_filename, vtk_dir, preserve_floodplain,
                  floodplain_elevation, do_inject_bathymetry):
    step = 5

    if do_inject_bathymetry:
        print('Step {}. Injecting bathymetry'.format(step))
        inject_bathymetry(mesh_file=out_filename)
        step += 1

    if preserve_floodplain:
        print('Step {}. Injecting flag to preserve floodplain'.format(step))
        inject_preserve_floodplain(mesh_file=out_filename,
                                   floodplain_elevation=floodplain_elevation)
        step += 1

    print('Step {}. Create vtk file for visualization'.format(step))
    extract_vtk(ignore_time=True, lonlat=True, dimension_list=['maxEdges='],
                variable_list=['allOnCells'], filename_pattern=out_filename,
                out_dir=vtk_dir)

    print("***********************************************")
    print("**    The global mesh file is {}   **".format(out_filename))
    print("***********************************************")
