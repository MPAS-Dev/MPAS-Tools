.. _ocean_coastal_tools:

.. |---| unicode:: U+2014  .. em dash, trimming surrounding whitespace
   :trim:

Coastal Tools
=============

The :py:mod:`mpas_tools.ocean.coastal_tools` module contains several functions
that aid in the construction of meshes refined around coastlines.

.. _coastal_tools_refine_mesh:

Refining a Mesh
---------------

The main driver of coastal tools is the function
:py:func:`mpas_tools.ocean.coastal_tools.coastal_refined_mesh()`.  This function
is called repeatedly with different dictionaries of ``params`` to add refinement
in different locations, starting from a background mesh (see
:ref:`coastal_tools_background_mesh`).

``coastal_tools`` provides the following dictionary of default parameters as
a starting point for ``params``:

.. code-block:: python

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
        # 'EC' (defaults to 60to30), 'QU' (uses dx_max_global),
        # 'RRS' (uses dx_max_global and dx_min_global)
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
        "nn_search": "flann",
        "plot_option": True

    }

A coastal mesh is defined by:

1. :ref:`coastal_tools_background_mesh`
2. :ref:`coastal_tools_extract`
3. :ref:`coastal_tools_distance`
4. :ref:`coastal_tools_cell_width`

Steps 2-4 can be repeated with several different regions to add more refinement.

The first time ``coastal_refined_mesh()`` is called, no cell width, latitude or
longitude coordinates are supplied.  In this case, a background mesh resolution
is defined through a call to
:py:func:`mpas_tools.ocean.coastal_tools.create_background_mesh()`.  The
following entries in ``params`` are used as arguments:

* ``'grd_box'`` - the bounds of the grid defining the mesh resolution
* ``'ddeg'`` - the resolution in longitude and latitude in degrees
* ``'mesh_type'`` - one of {``'QU'``, ``'EC'`` or ``'RRS'``} indicating the
  type of mesh
* ``'dx_min_global'`` - the resolution for a ``QU`` mesh, the minimum resolution
  for an ``RRS`` mesh, ignored for ``EC`` meshes.
* ``'dx_max_global'`` - the maximum resolution for an ``RRS`` mesh, ignored for
  ``QU`` and ``EC`` meshes.
* ``'plot_option'`` - whether to plot the results and save it to a PNG file.
* ``'plot_box'`` - If ``plot_option`` is ``True``, the bounds of the plot in
  longitude and latitude.

After the background field of cell widths has been created or using the
cell widths passed in as function arguments, ``coastal_refined_mesh()`` then
adds a region of refinement.

First, a set of coastline contours is extracted by calling
:py:func:`mpas_tools.ocean.coastal_tools.extract_coastlines()` with the
following values from ``params``:

* ``'nc_file'`` - A bathymetry dataset on a lon/lat grid in NetCDF format
* ``'nc_vars'`` - The names of the lon, lat and bathymetry variables in a list
* ``'region_box'`` - see :ref:`coastal_tools_regions`
* ``'z_contour'`` - A contour level from the bathymetry dataset to extract as a
  region boundary
* ``'n_longest'`` - The maximum number of contours to retain (sorted from
  longest to shortest).
* ``'point_list'`` - An optional list of points to add to the coastline
* ``'plot_option'`` - Whether to plot the extracted coastline.
* ``'plot_box'`` - If ``plot_option`` is ``True``, the bounds of the plot in
  longitude and latitude.

Next, the distance to the coastal contours is computed using
:py:func:`mpas_tools.ocean.coastal_tools.distance_to_coast()` with the
following values from ``params``:

* ``'origin'`` - A lon and lat point---no longer used in the code
* ``'nn_search'`` - Whether to use the ``'flann'`` or ``'kdtree'`` algorithm,
  with the ``'flann'`` strongly recommended.
* ``'smooth_coastline'`` - The number of neighboring cells along the coastline
  over which to average locations to smooth the coastline
* ``'plot_option'`` - Whether to plot the distance function.
* ``'plot_box'`` - If ``plot_option`` is ``True``, the bounds of the plot in
  longitude and latitude.

Finally, the distance function is used to blend the background and refined
regions using
:py:func:`mpas_tools.ocean.coastal_tools.compute_cell_width()` with the
following values from ``params``:

* ``'dx_min_coastal'`` - the resolution in meters of the refined region
* ``'trans_start'`` - the distance in meters from the coast at which the
  transition in resolution begins---the center of the transition is half a
  ``trans_width`` farther from the coastline
* ``'trans_width'`` - the distance in meters over which the transition occurs
* ``'restrict_box'`` - A region of made up of quadrilaterals to ``include`` and
  ``exclude`` that defines where resolution may be altered.  Outside of the
  ``restrict_box``, the resolution remains unchanged.  See
  :ref:`coastal_tools_regions`.
* ``'plot_option'`` - Whether to plot the cell widths and transition function.
* ``'plot_box'`` - If ``plot_option`` is ``True``, the bounds of the plot in
  longitude and latitude.

Here is an example of multiple calls to ``coastal_refined_mesh()`` in action,
taken from the
`hurricane/USDEQU120at30cr10rr2 <https://github.com/MPAS-Dev/MPAS-Model/blob/ocean/develop/testing_and_setup/compass/ocean/hurricane/USDEQU120at30cr10rr2/build_mesh/build_base_mesh.py>`_
test case from
`COMPASS <https://github.com/MPAS-Dev/MPAS-Model/tree/ocean/develop/testing_and_setup/compass>`_.
This workflow refines a background uniform mesh with 120-km resolution with
successively higher and higher resolution down to the Delaware Bay at 2-km
resolution.

.. code-block:: python

    import mpas_tools.ocean.coastal_tools as ct


    km = 1000.0

    params = ct.default_params

    print("****QU 120 background mesh and enhanced Atlantic (30km)****")
    params["mesh_type"] = "QU"
    params["dx_max_global"] = 120.0 * km
    params["region_box"] = ct.Atlantic
    params["restrict_box"] = ct.Atlantic_restrict
    params["plot_box"] = ct.Western_Atlantic
    params["dx_min_coastal"] = 30.0 * km
    params["trans_width"] = 5000.0 * km
    params["trans_start"] = 500.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(params)

    print("****Northeast refinement (10km)***")
    params["region_box"] = ct.Delaware_Bay
    params["plot_box"] = ct.Western_Atlantic
    params["dx_min_coastal"] = 10.0 * km
    params["trans_width"] = 600.0 * km
    params["trans_start"] = 400.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

    print("****Delaware regional refinement (5km)****")
    params["region_box"] = ct.Delaware_Region
    params["plot_box"] = ct.Delaware
    params["dx_min_coastal"] = 5.0 * km
    params["trans_width"] = 175.0 * km
    params["trans_start"] = 75.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

    print("****Delaware Bay high-resolution (2km)****")
    params["region_box"] = ct.Delaware_Bay
    params["plot_box"] = ct.Delaware
    params["restrict_box"] = ct.Delaware_restrict
    params["dx_min_coastal"] = 2.0 * km
    params["trans_width"] = 100.0 * km
    params["trans_start"] = 17.0 * km

    cell_width, lon, lat = ct.coastal_refined_mesh(
        params, cell_width, lon, lat)

.. _coastal_tools_background_mesh:

Creating a Background Mesh
--------------------------

A background mesh is typically created by calling
:py:func:`mpas_tools.ocean.coastal_tools.coastal_refined_mesh()` without
providing an input mesh but can also be created by calling
:py:func:`mpas_tools.ocean.coastal_tools.create_background_mesh()` directly.

The user must define the bounds of the grid in longitude and latitude (typically
-180 to 180 and -90 to 90, respectively) and the resolution in degrees.  The
mesh can be any of three types: {``'QU'``, ``'EC'`` or ``'RRS'``}.  For
Quasi-Uniform (QU) meshes, the resulting cell width will be a constant equal to
``dx_min_global``.  For Eddy-Closure (EC) meshes, the default parameters are
always used (see :ref:`ec_mesh`).  For Rossby-Radius Scaled (RRS) meshes,
``dx_min_global`` is the resolution at the poles while ``dx_max_global`` is the
resolution at the equator (see :ref:`rrs_mesh`).

.. _coastal_tools_extract:

Extracting Coastlines
---------------------

``coastal_tools`` extracts points along a coastline using
:py:func:`mpas_tools.ocean.coastal_tools.extract_coastlines()`.  The default
parameters are set up to use the
`earth_relief_15s.nc dataset <https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-ocean/bathymetry_database/SRTM15_plus_earth_relief_15s.nc>`_,
but any bathymetry data set on a lon/lat grid could be used as long as
``params['nc_file']`` is modified to point to the new name of the dataset and
``params['nc_vars']`` is set to the appropriate variable names.

By default, the coastline is extracted using the ``z = 0.0`` contour of the
bathyemtry but other values can be selected (e.g. to use distance from the
continental shelf break) by defining ``params['z_contour']``.
By default, only the 10 longest contours are retained to reduce computational
cost but more (or fewer) contours can be retained by setting
``params['n_longest']`` to another number.

Optionally, the results can be plotted withing the given "plot box" and saved
to a file.

.. _coastal_tools_distance:

Computing Distance to Coast
---------------------------

A key ingredient in defining resolution in coastal meshes is a field containing
the distance from each location in the field to the nearest point on the
coastline.  This distance field ``D`` is computed with
:py:func:`mpas_tools.ocean.coastal_tools.distance_to_coast()`
The user can optionally control the search algorithm used via
``params['nn_search']`` (though ``'flann'``, the default, is highly
recommended).  They can also decide to smooth the coastline as long as there is
a single coastline contour---with multiple contours, the current algorithm will
average the end of one contour with the start fo the next---by specifying an
integer number of neighbors as ``params['smooth_coastline']``.  The default is
no smoothing (``0`` neighbors).

.. _coastal_tools_cell_width:

Blending Cell Widths
--------------------

The final step in each iteration of coastal refinement is to blend the new,
refined resolution into the previous grid of cell widths with the function
:py:func:`mpas_tools.ocean.coastal_tools.compute_cell_width()`.

The most important parameters to set are ``params['dx_min_coastal']``, the
resolution in meters of the refined region; ``params['trans_start']``, the
distance from the coast in meters where the transition in resolution should
start; and ``params['trans_width']``, the width of the transition itself in
meters.

The resolution refinement can be confined to a region using
``params['restrict_box']`` to supply a region of made up of quadrilaterals to
``include`` and ``exclude`` from the restricted region. ``include`` boxes
specify regions where the distance-based coastal refinement function should be
used to update the desired cell width values. ``exclude`` boxes eliminate areas
within the ``include`` regions so they are not affected by the coastal
refinement function. This is useful for preventing resolution from appearing on
one side of a narrow piece of land. For example, coastal refinement in the Gulf
of Mexico should not appear on the other Pacific Ocean side of Central America.
The ``params['restrict_box']`` regions can be used to enforce this.

Here is an example of a restriction box for the region around Delaware, used in
the example above:

.. code-block:: python

    Delaware_restrict = {"include": [np.array([[-75.853, 39.732],
                                               [-74.939, 36.678],
                                               [-71.519, 40.156],
                                               [-75.153, 40.077]]),
                                     np.array([[-76.024, 37.188],
                                               [-75.214, 36.756],
                                               [-74.512, 37.925],
                                               [-75.274, 38.318]])],
                         "exclude": []}


.. _coastal_tools_regions:

Regions
-------

:ref:`coastal_tools_extract` requires a set of bounding regions to be defined.
These regions are made up of a list of quadrilaterals to ``include`` and another
list to ``exclude``.  The quadrilaterals are either bounding rectangles
(min lon, max lon, min lat, max lat) or lists of 4 (lon, lat) points numbered
counter clockwise.

``include`` boxes specify areas where coastline contours will be extracted from
the underlying bathymetry dataset. ``exclude`` boxes can be used to select
regions within ``include`` regions where coastline extraction should not be
done. For example, a large include box covering the U.S. East Coast may contain
small islands in the Atlantic, which may need to be ignored when placing coastal
refinement. In this case, ``exclude`` boxes can be specified to eliminate the
small coastlines.

An example of such a region is:

.. code-block:: python

    Greenland = {"include":[np.array([-81.5, -12.5, 60, 85])],
                 "exclude":[np.array([[-87.6, 58.7],
                                      [-51.9, 56.6],
                                      [-68.9, 75.5],
                                      [-107.0, 73.3]]),
                            np.array([[-119.0, 74.5],
                                      [-92.7, 75.9],
                                      [-52.84, 83.25],
                                      [-100.8, 84.0]]),
                            np.array([[-101.3, 68.5],
                                      [-82.4, 72.3],
                                      [-68.7, 81.24],
                                      [-117.29, 77.75]]),
                            np.array([-25.0, -10.0, 62.5, 67.5])]}

.. note::

    This example includes both bounding rectangles (e.g.
    ``np.array([-81.5, -12.5, 60, 85])``) and more general quadrilaterals (e.g.
    ``np.array([[-101.3, 68.5], [-82.4, 72.3],...``)

``coastal_tools`` defines 16 regions via dictionaries of this type.  The defined
regions are:

* Delaware_Bay
* Galveston_Bay
* Delaware_Region
* US_East_Coast
* US_Gulf_Coast
* Caribbean
* US_West_Coast
* Hawaii
* Alaska
* Bering_Sea_E
* Bering_Sea_W
* Aleutian_Islands_E
* Aleutian_Islands_W
* Greenland
* CONUS
* Continental_US

