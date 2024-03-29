#!/usr/bin/env python
"""
name: coastal_tools
authors: Steven Brus

last modified: 07/09/2018

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
from skimage import measure
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from scipy.io import savemat
import timeit
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import mpas_tools.mesh.creation.mesh_definition_tools as mdt
from mpas_tools.mesh.creation.util import lonlat2xyz

# Constants
km = 1000.0
deg2rad = np.pi / 180.0
rad2deg = 180.0 / np.pi

call_count = 0

##########################################################################
# Bounding box declarations (for coastal refinements)
##########################################################################

# ---------------
# Region boxes
# ---------------

# Bays and estuaries
Delaware_Bay = {"include": [np.array([-75.61903, -74.22, 37.8, 40.312747])],
                "exclude": []}
Galveston_Bay = {"include": [np.array([-95.45, -94.4, 29, 30])],
                 "exclude": []}

# Regions
Delaware_Region = {"include": [np.array([-77, -69.8, 35.5, 41])],
                   "exclude": []}

# Coastlines
US_East_Coast = {"include": [np.array([-81.7, -62.3, 25.1, 44.50])],  # East FL to Bay of Fundy
                 "exclude": [np.array([-66.0, -64.0, 31.5, 33.0]),    # Bermuda
                             np.array([-79.75, -70.0, 20.0, 28.5]),   # Bahamas
                             np.array([-65.15, -62.43, 43.0, 45.55]), # Gulf of St. Lawence
                             np.array([-66.65, -62.43, 43.0, 45.0])]} # ''

US_Gulf_Coast = {"include": [np.array([-98.0, -80.0, 24.0, 31.0]),    # West FL to NE Mexico
                             np.array([-98.5, -95.5, 20.0, 24.0]),    # East Mexico
                             np.array([-91.0, -86.0, 20.0, 22.0])     # Yucatan
                             ],
                 "exclude": []}

Caribbean = {"include": [np.array([-89.85, -59.73, 9.35, 20.86]),
                         ],
             "exclude": []}

US_West_Coast = {"include": [np.array([-127.0, -116.0, 32.5, 49.0]),  # California
                             np.array([-117.5, -108.0, 22.8, 32.5])   # Baja and West Mexico
                             ],
                 "exclude": [np.array([-116.5, -115.0, 32.8, 33.8]),  # Salton Sea
                             np.array([-120.5, -116.5, 35.5, 40.5]),  # Lake Tahoe, etc.
                             np.array([[-111.89, 21.24],              # Baja
                                       [-107.17, 22.48],
                                       [-113.94, 30.77],
                                       [-119.44, 33.09]])]}

Hawaii = {"include": [np.array([-161.0, -154.0, 18.0, 23.0])],
          "exclude": []}

Alaska = {"include": [np.array([-170.0, -141.0, 49.0, 63.0]),
                      np.array([-141.0, -129.5, 49.0, 61.0]),  # Southeast
                      np.array([-129.5, -121.0, 49.0, 55.0])   # Connects AK to CA
                      ],
          "exclude": [np.array([-171.0, -151.79, 49.54, 58.83])]}  # Aleutian Islands

Bering_Sea_E = {"include": [np.array([-180.0, -168.5, 56.00, 64.0])],
                "exclude": []}

Bering_Sea_W = {"include": [np.array([161.90, 180.0, 57.85, 66.63])],
                "exclude": []}

Aleutian_Islands_E = {"include": [np.array([-180.0, -151.79, 49.54, 58.83])],
                      "exclude": [np.array([-173.16, -164.37, 55.81, 57.55])]}

Aleutian_Islands_W = {"include": [np.array([164.9, 180.0, 50.77, 56.31])],
                      "exclude": [np.array([178.5, 179.5, 51.25, 51.75])]}

Greenland = {"include":[np.array([-81.5,-12.5,60,85])],
             "exclude":[np.array([[-87.6,58.7],
                                  [-51.9,56.6],
                                  [-68.9,75.5],
                                  [-107.0,73.3]]),
                        np.array([[-119.0,74.5],
                                  [-92.7,75.9],
                                  [-52.84,83.25],
                                  [-100.8,84.0]]),
                        np.array([[-101.3,68.5],
                                  [-82.4,72.3],
                                  [-68.7,81.24],
                                  [-117.29,77.75]]),
                        np.array([-25.0,-10.0,62.5,67.5])]}
Atlantic = {"include": [np.array([-78.5, -77.5, 23.5, 25.25])],  # Bahamas, use with large transition start to fill Atlantic
            "exclude": []}

# Combined coastlines
CONUS = {"include": [], "exclude": []}
CONUS["include"].extend(US_East_Coast["include"])
CONUS["include"].extend(US_Gulf_Coast["include"])
CONUS["include"].extend(US_West_Coast["include"])
CONUS["exclude"].extend(US_East_Coast["exclude"])
CONUS["exclude"].extend(US_West_Coast["exclude"])

Continental_US = {"include": [], "exclude": []}
Continental_US["include"].extend(CONUS["include"])
Continental_US["include"].extend(Alaska["include"])
Continental_US["exclude"].extend(CONUS["exclude"])

# ----------------
# Plotting boxes
# ----------------

Western_Atlantic = np.array([-98.186645, -59.832744, 7.791301, 45.942453])
Contiguous_US = np.array([-132.0, -59.832744, 7.791301, 51.0])
North_America = np.array([-175.0, -60.0, 7.5, 72.0])
Delaware = np.array([-77, -69.8, 35.5, 41])
Entire_Globe = np.array([-180, 180, -90, 90])

# -----------------
# Restrict Boxes
# -----------------

Empty = {"include": [],
         "exclude": []}

Delaware_restrict = {"include": [np.array([[-75.853, 39.732],
                                           [-74.939, 36.678],
                                           [-71.519, 40.156],
                                           [-75.153, 40.077]]),
                                 np.array([[-76.024, 37.188],
                                           [-75.214, 36.756],
                                           [-74.512, 37.925],
                                           [-75.274, 38.318]])],
                     "exclude": []}

Gulf_restrict = {"include": [np.array([[-85.04, 13.80],
                                       [-76.90, 16.60],
                                       [-86.24, 36.80],
                                       [-105.55, 22.63]])],
                 "exclude": []}

Caribbean_restrict = {"include": [np.array([[-76.39, 4.55],
                                            [-53.22, 4.29],
                                            [-53.22, 38.94],
                                            [-94.99, 18.47]])],
                      "exclude": [np.array([[-80.72, 1.66],
                                            [-73.70, 3.03],
                                            [-78.94, 9.33],
                                            [-84.98, 7.67]]),
                                  np.array([[-100.18, 13.76],
                                            [-82.93, 6.51],
                                            [-85.08, 13.74],
                                            [-95.86, 18.04]])]}

East_Coast_restrict = {"include": [],
                       "exclude": [np.array([[-72.0, 46.69],
                                             [-61.74, 45.48],
                                             [-44.07, 49.49],
                                             [-63.47, 53.76]])]}
Bering_Sea_restrict = {"include": [],
                       "exclude": [np.array([[143.46, 51.79],
                                             [155.65, 50.13],
                                             [166.23, 63.78],
                                             [148.63, 62.30]]),
                                   np.array([[154.36, 68.22],
                                             [-173.80, 65.94],
                                             [-161.81, 72.02],
                                             [163.64, 73.70]])]}

Atlantic_restrict = {"include": [np.array([[-86.39, 13.67],
                                           [-24.44, 21.32],
                                           [-50.09, 55.98],
                                           [-105.58, 23.61]]),
                                 np.array([[-76.39, 4.55],
                                           [-30.74, -2.66],
                                           [-30.83, 43.95],
                                           [-94.99, 18.47]])],
                     "exclude": [np.array([[-80.72, 1.66],
                                           [-73.70, 3.03],
                                           [-78.94, 9.33],
                                           [-84.98, 7.67]]),
                                 np.array([[-100.18, 13.76],
                                           [-82.93, 6.51],
                                           [-85.08, 13.74],
                                           [-95.86, 18.04]])]}

##########################################################################
# User-defined inputs
##########################################################################

default_params = {

    # Path to bathymetry data and name of file
    "nc_file": "./earth_relief_15s.nc",
    "nc_vars": ["lon","lat","z"],

    # Bounding box of coastal refinement region
    "region_box": Continental_US,
    "origin": np.array([-100, 40]),
    "restrict_box": Empty,

    # Coastline extraction parameters
    "z_contour": 0.0,
    "n_longest": 10,
    "smooth_coastline": 0,
    "point_list": None,

    # Global mesh parameters
    "grd_box": Entire_Globe,
    "ddeg": .1,
    # 'EC' (defaults to 60to30), 'QU' (uses dx_max_global), 'RRS' (uses dx_max_global and dx_min_global)
    "mesh_type": 'EC',
    "dx_max_global": 30.0 * km,
    "dx_min_global": 10.0 * km,

    # Coastal mesh parameters
    "dx_min_coastal": 10.0 * km,
    "trans_width": 600.0 * km,
    "trans_start": 400.0 * km,

    # Bounding box of plotting region
    "plot_box": North_America,

    # Options
    "nn_search": "kdtree",
    "plot_option": True

}

##########################################################################
# Functions
##########################################################################


def coastal_refined_mesh(params, cell_width=None, lon_grd=None, lat_grd=None):
    # {{{
    """
    Optionally create a background field of cell widths, then add a region of
    refined resolution to the cell widths.

    Parameters
    ----------
    params : dict
        A dictionary of parameters determining how the mesh is constructed.
        See ``mpas_tools.ocean.coastal_tools.default_params``.

    cell_width : ndarray, optional
        A 2D array of cell widths in meters.  If none is provided, one a base
        ``cell_width`` field constructed using parameter values from ``params``
        to call ``create_background_mesh``.

    lon_grd : ndarray, optional
        A 1D array of longitudes in degrees in the range from -180 to 180

    lat_grd : ndarray, optional
        A 1D array of latitudes in degrees in the range from -90 to 90

    Returns
    -------
    cell_width : ndarray
        A 2D array of cell widths in meters.

    lon_grd : ndarray
        A 1D array of longitudes in degrees in the range from -180 to 180

    lat_grd : ndarray
        A 1D array of latitudes in degrees in the range from -90 to 90
    """

    coastal_refined_mesh.counter += 1
    call_count = coastal_refined_mesh.counter

    if cell_width is None:
        # Create the background cell width array
        lon_grd, lat_grd, cell_width = create_background_mesh(
            params["grd_box"],
            params["ddeg"],
            params["mesh_type"],
            params["dx_min_global"],
            params["dx_max_global"],
            params["plot_option"],
            params["plot_box"],
            call_count)

    # Get coastlines from bathy/topo  data set
    coastlines = extract_coastlines(
        params["nc_file"],
        params["nc_vars"],
        params["region_box"],
        params["z_contour"],
        params["n_longest"],
        params["point_list"],
        params["plot_option"],
        params["plot_box"],
        call_count)

    # Compute distance from background grid points to coastline
    D = distance_to_coast(
        coastlines,
        lon_grd,
        lat_grd,
        params["nn_search"],
        params["smooth_coastline"],
        params["plot_option"],
        params["plot_box"],
        call_count)

    # Blend coastline and background resolution, save cell_width array as .mat file
    cell_width = compute_cell_width(
        D,
        cell_width,
        lon_grd,
        lat_grd,
        params["dx_min_coastal"],
        params["trans_start"],
        params["trans_width"],
        params["restrict_box"],
        params["plot_option"],
        params["plot_box"],
        coastlines,
        call_count)

    # Save matfile
    # save_matfile(cell_width/km,lon_grd,lat_grd)
    print("")

    return (cell_width, lon_grd, lat_grd)


    # }}}
coastal_refined_mesh.counter = 0

##############################################################


def create_background_mesh(grd_box, ddeg, mesh_type, dx_min, dx_max,  # {{{
                           plot_option=False, plot_box=[], call=None):
    """
    Create a background field of cell widths

    Parameters
    ----------
    grd_box : list of float
        A list of 4 floats defining the bounds (min lon, max lon, min lat, max
        lat) of the grid

    ddeg : float
        The resolution of the mesh in degrees

    mesh_type : {'QU', 'EC', 'RRS'}
        The type of mesh: quasi-uniform (QU), Eddy-closure (EC) or Rossby-radius
        scaling (RRS)

    dx_min : float
        The resolution in meters of a QU mesh or the minimum resolution of of
        an RRS mesh. This parameter is ignored for EC meshes and the default
        function arguments to ``EC_CellWidthVsLat()`` are used instead.

    dx_max : float
        The maximum resolution in meters of of an RRS mesh. This parameter is
        ignored for QU meshes and EC meshes.  For EC meshes, the default
        function arguments are used instead.

    plot_option : bool, optional
        Whether to plot the resulting cell width and save it to files
        named ``bckgrnd_grid_cell_width_vs_lat###.png`` and
        ``bckgnd_grid_cell_width###.png``, where ``###`` is given by
        ``call`` and is meant to indicate how many times this function has been
        called during mesh creation.

    plot_box : list of float, optional
        The extent of the plot if ``plot_option=True``

    call : int, optional
        The number of times the function has been called, used to give the
        plot a unique name.


    Returns
    -------
    cell_width : ndarray
        A 2D array of cell widths in meters.

    lon_grd : ndarray
        A 1D array of longitudes in degrees in the range from -180 to 180

    lat_grd : ndarray
        A 1D array of latitudes in degrees in the range from -90 to 90
    """

    print("Create background mesh")
    print("------------------------")

    # Create cell width background grid
    ny_grd = int((grd_box[3]-grd_box[2])/ddeg) + 1
    nx_grd = int((grd_box[1]-grd_box[0])/ddeg) + 1
    lat_grd = grd_box[2] + ddeg*np.arange(ny_grd)
    lon_grd = grd_box[0] + ddeg*np.arange(nx_grd)
    print("   Background grid dimensions:", ny_grd, nx_grd)

    # Assign background grid cell width values
    if mesh_type == 'QU':
        cell_width_lat = dx_max / km * np.ones(lat_grd.size)
    elif mesh_type == 'EC':
        cell_width_lat = mdt.EC_CellWidthVsLat(lat_grd)
    elif mesh_type == 'RRS':
        cell_width_lat = mdt.RRS_CellWidthVsLat(lat_grd, dx_max / km, dx_min / km)
    else:
        raise ValueError('Unknown mesh_type {}'.format(mesh_type))
    cell_width = np.tile(cell_width_lat, (nx_grd, 1)).T * km

    # Plot background cell width
    if plot_option:
        print("   Plotting background cell width")

        plt.figure()
        plt.plot(lat_grd, cell_width_lat)
        plt.savefig('bckgrnd_grid_cell_width_vs_lat' + str(call) + '.png')

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        plt.contourf(lon_grd, lat_grd, cell_width, transform=ccrs.PlateCarree())
        plot_coarse_coast(ax, plot_box)
        plt.colorbar()
        plt.savefig(
            'bckgnd_grid_cell_width' +
            str(call) +
            '.png',
            bbox_inches='tight')
        plt.close()
        print("   Done")

    return (lon_grd, lat_grd, cell_width)  # }}}

##############################################################


def extract_coastlines(nc_file, nc_vars, region_box, z_contour=0, n_longest=10,
                       point_list=None, plot_option=False, plot_box=None,
                       call=None):   # {{{

    """
    Extracts a set of coastline contours

    Parameters
    ----------
    nc_file : str
        A bathymetry dataset on a lon/lat grid in NetCDF format

    nc_vars : list of str
        The names of the longitude (nc_vars[0]), latitude (nc_vars[1]) and
        bathymetry (nc_vars[2]) variables.

    region_box : dict of list of ndarrays
        A region made up of a list of quadrilaterals to ``include`` and another
        list to ``exclude``.  The quadrilaterals are either bounding rectangles
        (min lon, max lon, min lat, max lat) or lists of 4 (lon, lat) points.

    z_contour : float, optional
        The isocontour of the bathymetry dataset to extract

    n_longest : int, optional
        The maximum number of contours to keep, after sorting from the longest
        to the shortest

    point_list : ndarray, optional
        A list of points to add to the coastline

    plot_option : bool, optional
        Whether to plot the resulting coastline points and the plot to a file
        named ``bathy_coastlines###.png``, where ``###`` is given by
        ``call`` and is meant to indicate how many times this function has been
        called during mesh creation.

    plot_box : list of float, optional
        The extent of the plot if ``plot_option=True``

    call : int, optional
        The number of times the function has been called, used to give the
        plot a unique name.

    Returns
    -------
    coastlines : ndarray
        An n x 2 array of (longitude, latitude) points along the coastline
        contours
    """

    print("Extract coastlines")
    print("------------------")

    # Open NetCDF file and read cooordintes
    nc_fid = Dataset(nc_file, "r")
    lon = nc_fid.variables[nc_vars[0]][:]
    lat = nc_fid.variables[nc_vars[1]][:]
    bathy_data = nc_fid.variables[nc_vars[2]]

    # Get coastlines for refined region
    coastline_list = []
    for i,box in enumerate(region_box["include"]):

        # Find coordinates and data inside bounding box
        xb,rect= get_convex_hull_coordinates(box)
        lon_region, lat_region, z_region = get_data_inside_box(lon, lat, bathy_data, xb)
        z_data = np.zeros(z_region.shape)
        z_data.fill(np.nan)
        idx = get_indices_inside_quad(lon_region,lat_region,box)
        z_data[idx] = z_region[idx]
        print("   Regional bathymetry data shape:", z_region.shape)

        # Find coastline contours
        print("   Extracting coastline "+str(i+1)+"/"+str(len(region_box["include"])))
        contours = measure.find_contours(z_data, z_contour)

        # Keep only n_longest coastlines and those not within exclude areas
        contours.sort(key=len, reverse=True)
        for c in contours[:n_longest]:
            # Convert from pixel to lon,lat
            c[:, 0] = (xb[3] - xb[2]) / float(len(lat_region)) * c[:, 0] + xb[2]
            c[:, 1] = (xb[1] - xb[0]) / float(len(lon_region)) * c[:, 1] + xb[0]
            c = np.fliplr(c)

            exclude = False
            for area in region_box["exclude"]:
                # Determine coastline coordinates in exclude area
                idx = get_indices_inside_quad(
                    c[:, 0], c[:, 1], area, grid=False)

                # Exlude coastlines that are entirely contained in exclude area
                if idx.size == c.shape[0]:
                    exclude = True
                    break
                elif idx.size != 0:
                    c = np.delete(c, idx, axis=0)

            # Keep coastlines not entirely contained in exclude areas
            if not exclude:
                cpad = np.vstack((c, [np.nan, np.nan]))
                coastline_list.append(cpad)

        print("   Done")

    # Add in user-specified points
    if point_list:
        for i,points in enumerate(point_list):
            cpad = np.vstack((points, [np.nan, np.nan]))
            coastline_list.append(cpad)

    # Combine coastlines
    coastlines = np.concatenate(coastline_list)

    if plot_option:
        print("   Plotting coastlines")

        # Find coordinates and data inside plotting box
        lon_plot, lat_plot, z_plot = get_data_inside_box(
            lon, lat, bathy_data, plot_box)

        # Plot bathymetry data, coastlines and region boxes

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        levels = np.linspace(np.amin(z_plot), np.amax(z_plot), 100)
        ds = 100                               # Downsample
        dsx = np.arange(0, lon_plot.size, ds)  # bathy data
        dsy = np.arange(0, lat_plot.size, ds)  # to speed up
        dsxy = np.ix_(dsy, dsx)                # plotting
        plt.contourf(lon_plot[dsx], lat_plot[dsy], z_plot[dsxy], levels=levels,
                     transform=ccrs.PlateCarree())
        plot_coarse_coast(ax, plot_box)
        plt.plot(coastlines[:, 0], coastlines[:, 1], color='white')
        for box in region_box["include"]:
            plot_region_box(box, 'b')
        for box in region_box["exclude"]:
            plot_region_box(box, 'r')
        plt.colorbar()
        plt.axis('equal')
        plt.savefig(
            'bathy_coastlines' +
            str(call) +
            '.png',
            bbox_inches='tight')
        plt.close()
        print("   Done")

    return coastlines  # }}}

##############################################################


def distance_to_coast(coastlines, lon_grd, lat_grd, nn_search='kdtree',
                      smooth_window=0, plot_option=False, plot_box=[],
                      call=None, workers=-1):
    # {{{
    """
    Extracts a set of coastline contours

    Parameters
    ----------
    coastlines : ndarray
        An n x 2 array of (longitude, latitude) points along the coastline
        contours returned from ``extract_coastlines()``

    lon_grd : ndarray
        A 1D array of longitudes in degrees in the range from -180 to 180

    lat_grd : ndarray
        A 1D array of latitudes in degrees in the range from -90 to 90

    nn_search : {'kdtree'}, optional
        The algorithm to use for the nearest neightbor search.

    smooth_window : int, optional
        The number of adjacent coastline points to average together to smooth
        the coastal contours.  Use ``0`` to indicate no smoothing.

    plot_option : bool, optional
        Whether to plot the resulting coastline points and the plot to a file
        named ``bathy_coastlines###.png``, where ``###`` is given by
        ``call`` and is meant to indicate how many times this function has been
        called during mesh creation.

    plot_box : list of float, optional
        The extent of the plot if ``plot_option=True``

    call : int, optional
        The number of times the function has been called, used to give the
        plot a unique name.

    workers : int, optional
        The number of threads used for finding nearest neighbors.  The default
        is all available threads (``workers=-1``)

    Returns
    -------
    D : ndarray
        A len(lat_grd) x len(lon_grd) array of distances in meters on the
        lon/lat grid to the closest point in the (smoothed) coastline contour.
    """

    if nn_search != 'kdtree':
        raise ValueError(f'nn_search method {nn_search} not available.')

    print("Distance to coast")
    print("-----------------")

    # Remove Nan values separating coastlines
    coast_pts = coastlines[np.isfinite(coastlines).all(axis=1)]

    # Smooth coast points if necessary
    if smooth_window > 1:
        coast_pts[:, 0], coast_pts[:, 1] = smooth_coastline(
            coast_pts[:, 0], coast_pts[:, 1], smooth_window)

    # Convert coastline points to x,y,z and create kd-tree
    npts = coast_pts.shape[0]
    coast_pts_xyz = np.zeros((npts,3))
    coast_pts_xyz[:, 0], coast_pts_xyz[:, 1], coast_pts_xyz[:, 2] = lonlat2xyz(coast_pts[:, 0], coast_pts[:, 1])
    tree = KDTree(coast_pts_xyz)

    # Convert  backgound grid coordinates to x,y,z and put in a nx_grd x 3 array for kd-tree query
    Lon_grd, Lat_grd = np.meshgrid(lon_grd, lat_grd)
    X_grd, Y_grd, Z_grd = lonlat2xyz(Lon_grd,Lat_grd)
    pts = np.vstack([X_grd.ravel(), Y_grd.ravel(), Z_grd.ravel()]).T

    # Find distances of background grid coordinates to the coast
    print("   Finding distance")
    start = timeit.default_timer()
    d, idx = tree.query(pts, workers=workers)
    end = timeit.default_timer()
    print("   Done")
    print("   " + str(end - start) + " seconds")

    # Make distance array that corresponds with cell_width array
    D = np.reshape(d, Lon_grd.shape)

    if plot_option:
        print("   Plotting distance to coast")

        # Find coordinates and data inside plotting box
        lon_plot, lat_plot, D_plot = get_data_inside_box(
            lon_grd, lat_grd, D, plot_box)

        # Plot distance to coast
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        D_plot = D_plot / km
        levels = np.linspace(np.amin(D_plot), np.amax(D_plot), 10)
        plt.contourf(lon_plot, lat_plot, D_plot, levels=levels,
                     transform=ccrs.PlateCarree())
        plot_coarse_coast(ax, plot_box)
        plt.plot(coastlines[:, 0], coastlines[:, 1], color='white')
        plt.grid(
            xdata=lon_plot,
            ydata=lat_plot,
            c='k',
            ls='-',
            lw=0.1,
            alpha=0.5)
        plt.colorbar()
        plt.axis('equal')
        plt.savefig('distance' + str(call) + '.png', bbox_inches='tight')
        plt.close()
        print("   Done")

    return D  # }}}

##############################################################


def compute_cell_width(D, cell_width, lon, lat, dx_min, trans_start,
                       trans_width, restrict_box, plot_option=False,
                       plot_box=[], coastlines=[], call=None):  # {{{
    """
    Blend cell widths from the input field with the new resolution in the
    refined region determined by the distance to the coastline contour.

    Parameters
    ----------
    D : ndarray
        A len(lat) x len(lon) array of distances in meters on the lon/lat grid
        to the closest point in the (smoothed) coastline contour returned from
        ``distance_to_coast()``

    cell_width : ndarray
        A len(lat) x len(lon) array of cell widths in meters

    lon : ndarray
        A 1D array of longitudes in degrees in the range from -180 to 180

    lat : ndarray
        A 1D array of latitudes in degrees in the range from -90 to 90

    dx_min : float
        The resolution in meters of the new refined region.

    trans_start : float
        The approximate value of ``D`` in meters at which the transition in
        resolution should start

    trans_width : float
        The approximate width in meters over which the transition in resolution
        should take place

    restrict_box : dict of lists of ndarrays
        A region of made up of quadrilaterals to ``include`` and ``exclude``
        that defines where resolution may be altered.  Outside of the
        ``restrict_box``, the resolution remains unchanged.

    plot_option : bool, optional
        Whether to plot the resulting coastline points and the plot to files
        named ``cell_width###.png`` and ``trans_func###.png```, where ``###``
        is given by ``call`` and is meant to indicate how many times this
        function has been called during mesh creation.

    plot_box : list of float, optional
        The extent of the plot if ``plot_option=True``

    coastlines : ndarray
        An n x 2 array of (longitude, latitude) points along the coastline
        contours returned from ``extract_coastlines()`` used in plotting if
        ``plot_option=True``

    call : int, optional
        The number of times the function has been called, used to give the
        plot a unique name.

    Returns
    -------
    cell_width : ndarray
        A len(lat) x len(lon) array of the new cell widths in meters
    """

    print("Compute cell width")
    print("------------------")

    # Compute cell width based on distance
    print("   Computing cell width")
    backgnd_weight = .5 * \
        (np.tanh((D - trans_start - .5 * trans_width) / (.2 * trans_width)) + 1.0)
    dist_weight = 1.0 - backgnd_weight
    ## Use later for depth and slope dependent resolution
    ##hres = np.maximum(dx_min*bathy_grd/20,dx_min)
    ##hres = np.minimum(hres,dx_max)
    #hw = np.zeros(Lon_grd.shape) + dx_max
    #hw[ocn_idx] = np.sqrt(9.81*bathy_grd[ocn_idx])*12.42*3600/25
    #hs = .20*1/dbathy_grd
    #h = np.fmin(hs,hw)
    #h = np.fmin(h,dx_max)
    #h = np.fmax(dx_min,h)

    cell_width_old = np.copy(cell_width)

    # Apply cell width function
    if len(restrict_box["include"]) > 0:
        # Only apply inside include regions
        for box in restrict_box["include"]:
            idx = get_indices_inside_quad(lon, lat, box)
            cell_width[idx] = (dx_min*dist_weight[idx] +
                               np.multiply(cell_width_old[idx], backgnd_weight[idx]))
    else:
        # Apply everywhere
        cell_width = (dx_min*dist_weight +
                      np.multiply(cell_width_old, backgnd_weight))

    # Don't apply cell width function in exclude regions (revert to previous values)
    if len(restrict_box["exclude"]) > 0:
        for box in restrict_box["exclude"]:
            idx = get_indices_inside_quad(lon, lat, box)
            cell_width[idx] = cell_width_old[idx]
    print("   Done")

    if plot_option:
        print("   Plotting cell width")

        # Find coordinates and data inside plotting box
        lon_plot, lat_plot, cell_width_plot = get_data_inside_box(
            lon, lat, cell_width / km, plot_box)

        # Plot cell width
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        levels = np.linspace(np.amin(cell_width_plot),
                             np.amax(cell_width_plot), 100)
        plt.contourf(lon_plot, lat_plot, cell_width_plot, levels=levels,
                     transform=ccrs.PlateCarree())
        plot_coarse_coast(ax, plot_box)
        plt.plot(coastlines[:, 0], coastlines[:, 1], color='white')
        if restrict_box:
            for box in restrict_box["include"]:
                plot_region_box(box, 'b')
            for box in restrict_box["exclude"]:
                plot_region_box(box, 'r')
        plt.colorbar()
        plt.axis('equal')
        plt.savefig('cell_width' + str(call) + '.png', bbox_inches='tight')
        plt.close()

        # Plot cell width transistion functions
        ts = trans_start/km
        tw = trans_width/km
        d = np.linspace(0,2*(ts+tw),1000)
        bw = .5*(np.tanh((d-ts-.5*tw)/(.2*tw))+1)
        dw = 1-bw
        plt.figure()
        plt.plot(d,bw)
        plt.plot(d,dw)
        plt.legend(('background','coastal region'))
        plt.plot([ts,ts],[0.0,1.0],'k-')
        plt.plot([ts+tw,ts+tw],[0.0,1.0],'k-')
        plt.tight_layout()
        plt.xlabel('distance (km)')
        plt.ylabel('weight')
        plt.savefig('trans_func'+str(call)+'.png',bbox_inches='tight')
        plt.close()
        print("   Done")

    return cell_width  # }}}

##############################################################


def save_matfile(cell_width, lon, lat):

    savemat('cellWidthVsLatLon.mat',
            mdict={'cellWidth': cell_width,
                   'lon': lon,
                   'lat': lat})

##############################################################


def CPP_projection(lon, lat, origin):

    R = 6378206.4
    origin = origin * deg2rad
    x = R * (lon * deg2rad - origin[0]) * np.cos(origin[1])
    y = R * lat * deg2rad

    return x, y

##############################################################


def smooth_coastline(x, y, window):
    xs = np.copy(x)
    ys = np.copy(y)
    offset = (window-1)/2
    for pt in range(offset-1, x.size-offset):
        xs[pt] = np.mean(x[pt-offset:pt+offset])
        ys[pt] = np.mean(y[pt-offset:pt+offset])
    return xs, ys
##############################################################


def get_data_inside_box(lon, lat, data, box, idx=False):

    # Find indicies of coordinates inside bounding box
    lon_idx, = np.where((lon >= box[0]) & (lon <= box[1]))
    lat_idx, = np.where((lat >= box[2]) & (lat <= box[3]))

    # Get indicies inside bounding box
    lon_region = lon[lon_idx]
    lat_region = lat[lat_idx]
    latlon_idx = np.ix_(lat_idx, lon_idx)

    # Return data inside bounding box
    if idx == False:

        try:     # Numpy indexing
            z_region = data[latlon_idx]
        except:  # NetCDF indexing
            z_region = data[lat_idx, lon_idx]

        return (lon_region, lat_region, z_region)

    # Return indicies of data inside bounding box
    elif idx == True:

        return latlon_idx

##############################################################


def get_indices_inside_quad(lon, lat, box, grid=True):

    wrap = flag_wrap(box)

    lon_adj = np.copy(lon)
    if wrap:
        idx = np.where((lon_adj >= -180.0) & (lon_adj <= -90.0))
        lon_adj[idx] = lon_adj[idx] + 360.0

    if grid:
        # Create vectors of all coordinates
        Lon_grd, Lat_grd = np.meshgrid(lon_adj, lat)
        X = Lon_grd.ravel()
        Y = Lat_grd.ravel()
    else:
        X = lon_adj
        Y = lat

    xb,rect = get_convex_hull_coordinates(box)

    # Find indices of coordinates in convex hull of quad region
    idxx = np.where((X >= xb[0]) & (X <= xb[1]))
    idxy = np.where((Y >= xb[2]) & (Y <= xb[3]))
    idx_ch = np.intersect1d(idxx, idxy)
    idx = np.copy(idx_ch)

    if rect == True:
        if grid == True:
            idx = np.unravel_index(idx, Lon_grd.shape)
        return idx

    # Initialize the local coordinate vectors to be outside unit square
    R = 0.0 * X - 10.0
    S = 0.0 * Y - 10.0

  # Initialize the coordinate vectors for points inside convex hull of quad region
    r = 0.0 * R[idx]
    s = 0.0 * S[idx]
    x = X[idx]
    y = Y[idx]

    # Map all coordinates in convex hull of quad region to unit square
    #   by solving inverse transformaion with Newton's method
    tol = 1e-8
    maxit = 25
    for it in range(0, maxit):

        # Compute shape fuctions
        phi1 = np.multiply((1.0 - r), (1.0 - s))
        phi2 = np.multiply((1.0 + r), (1.0 - s))
        phi3 = np.multiply((1.0 + r), (1.0 + s))
        phi4 = np.multiply((1.0 - r), (1.0 + s))

        # Compute functions that are being solved
        f1 = .25 * (phi1 * box[0, 0] + phi2 * box[1, 0] +
                    phi3 * box[2, 0] + phi4 * box[3, 0]) - x
        f2 = .25 * (phi1 * box[0, 1] + phi2 * box[1, 1] +
                    phi3 * box[2, 1] + phi4 * box[3, 1]) - y

        # Compute Jacobian
        df1ds = .25 * ((r - 1.0) * box[0, 0] - (1.0 + r) * box[1, 0] + (1.0 + r) * box[2, 0] + (1.0 - r) * box[3, 0])
        df1dr = .25 * ((s - 1.0) * box[0, 0] + (1.0 - s) * box[1, 0] + (1.0 + s) * box[2, 0] - (1.0 + s) * box[3, 0])
        df2ds = .25 * ((r - 1.0) * box[0, 1] - (1.0 + r) * box[1, 1] + (1.0 + r) * box[2, 1] + (1.0 - r) * box[3, 1])
        df2dr = .25 * ((s - 1.0) * box[0, 1] + (1.0 - s) * box[1, 1] + (1.0 + s) * box[2, 1] - (1.0 + s) * box[3, 1])

        # Inverse of 2x2 matrix
        det_recip = np.multiply(df1dr, df2ds) - np.multiply(df2dr, df1ds)
        det_recip = 1.0 / det_recip
        dr = np.multiply(det_recip,  np.multiply(df2ds, f1) - np.multiply(df1ds, f2))
        ds = np.multiply(det_recip, -np.multiply(df2dr, f1) + np.multiply(df1dr, f2))

        # Apply Newton's method
        rnew = r - dr
        snew = s - ds

        # Find converged values
        err = R[idx] - rnew
        idxr = np.where(np.absolute(err) < tol)
        err = S[idx] - snew
        idxs = np.where(np.absolute(err) < tol)
        idx_conv = np.intersect1d(idxr, idxs)

        # Update solution
        R[idx] = rnew
        S[idx] = snew

        # Find indicies of unconverged values
        idx = np.delete(idx, idx_conv)
        # print("Iteration: ", it, "unconverged values: ", idx.size)

        # Terminate once all values are converged
        if idx.size == 0:
            break

        # Initialize to unconverged values for next iteration
        r = R[idx]
        s = S[idx]
        x = X[idx]
        y = Y[idx]

    # Find any remaining unconverged values
    if grid == True:
        idx_nc = np.unravel_index(idx, Lon_grd.shape)
    else:
        idx_nc = np.copy(idx)

    # Find indicies of coordinates inside quad region
    lon_idx, = np.where((R >= -1.0) & (R <= 1.0))
    lat_idx, = np.where((S >= -1.0) & (S <= 1.0))
    idx = np.intersect1d(lon_idx, lat_idx)
    if grid == True:
        idx = np.unravel_index(idx, Lon_grd.shape)

    ## Plot values inside quad region
    # plt.figure()
    # plt.plot(X,Y,'.')
    # if grid == True:
    #  plt.plot(Lon_grd[idx],Lat_grd[idx],'.')
    #  plt.plot(Lon_grd[idx_nc],Lat_grd[idx_nc],'.')
    # else:
    #  plt.plot(lon[idx],lat[idx],'.')
    #  plt.plot(lon[idx_nc],lat[idx_nc],'.')
    # plt.plot(box[:,0],box[:,1],'o')
    # plt.savefig("restrict_box.png")

    return idx

##############################################################

def get_convex_hull_coordinates(box):

    wrap = flag_wrap(box)

    xb = np.zeros(4)
    if box.size == 4:
        if wrap:
            for i in range(2):
                if box[i] >= -180.0 and box[i] <= -90.0:
                    box[i] = box[i] + 360.0

        xb[0] = box[0]
        xb[1] = box[1]
        xb[2] = box[2]
        xb[3] = box[3]
        rect = True
    else:
        if wrap:
            for i in range(4):
                if box[i, 0] >= -180.0 and box[i, 0] <= -90.0:
                    box[i, 0] = box[i, 0] + 360.0

        xb[0] = np.amin(box[:, 0])
        xb[1] = np.amax(box[:, 0])
        xb[2] = np.amin(box[:, 1])
        xb[3] = np.amax(box[:, 1])
        rect = False

    return xb,rect

##############################################################

def flag_wrap(box):

    wrap = False
    if box.size == 4:
      if box[0] > 0.0 and box[1] < 0.0:
        wrap = True
    else:
      if np.any(box[:,0] > 0.0) and np.any(box[:,0] < 0.0):
        wrap = True

    return wrap

##############################################################

def plot_coarse_coast(ax, plot_box):

    ax.set_extent(plot_box, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)

##############################################################


def plot_region_box(box, color):
    ls = color + '-'
    if box.size == 4:
        plt.plot([box[0], box[1]], [box[2], box[2]], ls)
        plt.plot([box[1], box[1]], [box[2], box[3]], ls)
        plt.plot([box[1], box[0]], [box[3], box[3]], ls)
        plt.plot([box[0], box[0]], [box[3], box[2]], ls)
    else:
        plt.plot([box[0, 0], box[1, 0]], [box[0, 1], box[1, 1]], ls)
        plt.plot([box[1, 0], box[2, 0]], [box[1, 1], box[2, 1]], ls)
        plt.plot([box[2, 0], box[3, 0]], [box[2, 1], box[3, 1]], ls)
        plt.plot([box[3, 0], box[0, 0]], [box[3, 1], box[0, 1]], ls)


##########################################################################
#  Incorporate later for depth and slope dependent resolution
##########################################################################


##
## Interpolate bathymetry onto background grid
#Lon_grd = Lon_grd*deg2rad
#Lat_grd = Lat_grd*deg2rad
#bathy = inject_bathymetry.interpolate_bathymetry(data_path+nc_file,Lon_grd.ravel(),Lat_grd.ravel())
#bathy_grd = -1.0*np.reshape(bathy,(ny_grd,nx_grd))
#ocn_idx = np.where(bathy_grd > 0)
#
# if plot_option:
#  plt.figure()
#  levels = np.linspace(0,11000,100)
#  plt.contourf(lon_grd,lat_grd,bathy_grd,levels=levels)
#  plot_coarse_coast(plot_box)
#  plt.colorbar()
#  plt.axis('equal')
#  plt.savefig('bckgnd_grid_bathy.png',bbox_inches='tight')

## Interpolate bathymetry gradient onto background grid
#dbathy = inject_bathymetry.interpolate_bathymetry(data_path+nc_file,Lon_grd.ravel(),Lat_grd.ravel(),grad=True)
#dbathy = np.reshape(dbathy,(ny_grd,nx_grd))
#dbathy_grd = np.zeros((ny_grd,nx_grd))
#dbathy_grd[ocn_idx] = dbathy[ocn_idx]
#
# if plot_option:
#  plt.figure()
#  plt.contourf(lon_grd,lat_grd,1/dbathy_grd)
#  plot_coarse_coast(plot_box)
#  plt.colorbar()
#  plt.axis('equal')
#  plt.savefig('bckgnd_grid_bathy_grad.png',bbox_inches='tight')
