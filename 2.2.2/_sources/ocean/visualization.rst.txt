.. _ocean_visualization:

*************
Visualization
*************

.. _ocean_viz_transects:

Plotting Ocean Transects
========================

.. image:: ../images/south_atlantic_temperature_transect.png
   :width: 500 px
   :align: center

The function :py:func:`mpas_tools.ocean.viz.transect.plot_transect()` is
used to plot transects of various MPAS-Ocean variables. This function provide
high-level plotting capabilities for visualizing transects of MPAS-Ocean
data. It can create transect plots useful for visualizing test cases or
analyzing global simulation output and comparing with observations.

The function :py:func:`mpas_tools.ocean.viz.transect.plot_feature_transects()`
and the associated ``plot_ocean_transects`` command-line tool can be used to
plot transects of various MPAS-Ocean variables.  The arguments to the
command-line tool are:

.. code-block:: none

    $ plot_ocean_transects --help
    usage: plot_ocean_transects [-h] -g GEOJSON_FILENAME [-m MESH_FILENAME] -f
                                FILENAME [-v VARIABLE_LIST [VARIABLE_LIST ...]]
                                [-c COLORMAP] [--flip]

    options:
      -h, --help            show this help message and exit
      -g GEOJSON_FILENAME, --geojson GEOJSON_FILENAME
                            A geojson file with transects to plot
      -m MESH_FILENAME, --mesh MESH_FILENAME
                            An MPAS-Ocean mesh file.  If not specified, the MPAS-Ocean data file must contain the mesh.
      -f FILENAME, --file FILENAME
                            An MPAS-Ocean data file
      -v VARIABLE_LIST [VARIABLE_LIST ...], --variable_list VARIABLE_LIST [VARIABLE_LIST ...]
                            List of variables to plot.  All variables on cells in the data file is the default.
      -c COLORMAP, --colormap COLORMAP
                            A colormap to use for the plots, default depends on the field name.
      --flip                Flip the x axis for all transects
      --write_netcdf        Whether to write a NetCDF file for the transect in addition to the image
      --method METHOD       The type of interpolation to use in plots. Options are "flat" and "bilinear"


See `transects <https://github.com/MPAS-Dev/geometric_features/tree/main/geometric_data/ocean/transect>`_
from ``geometric_features`` for a examples of what a geojson transect might
look like:

.. code-block:: json

    {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "name": "Drake Passage",
                    "object": "transect",
                    "component": "ocean",
                    "author": "Mark Petersen, Xylar Asay-Davis, Milena Veneziani",
                },
                "geometry": {
                    "type": "LineString",
                    "coordinates": [
                        [
                            -63.02,
                            -65.46
                        ],
                        [
                            -63.81,
                            -63.8
                        ],
                        [
                            -64.42,
                            -62.02
                        ],
                        [
                            -65.04,
                            -60.25
                        ],
                        [
                            -65.74,
                            -58.28
                        ],
                        [
                            -66.37,
                            -56.39
                        ],
                        [
                            -67.02,
                            -54.44
                        ]
                    ]
                }
            }
        ]
    }

Add more features to the ``features`` list to plot multiple transects at the
same time.

The MPAS-Ocean mesh file must including not just the horizontal mesh variables
but also the vertical mesh variables (``minLevelCell``, ``maxLevelCell``,
``layerThickness``, etc.)

If you don't specify the list of variables to plot, all variables with
dimensions ``nCells`` and ``nVertLevels`` will be plotted.

One way of customizing these visualizaitons is to make your own copy of
`transects.py <https://github.com/MPAS-Dev/MPAS-Tools/blob/master/conda_package/mpas_tools/ocean/viz/transects.py>`_
and change ``_plot_transect()`` to suite your needs, (changing figure size, dpi,
colorbar, etc.)

.. _ocean_viz_transects_interp:

Ocean Transect Interpolation
============================

The ``mpas_tools.ocean.viz.transect.vert`` module provides functions for
interpolating MPAS-Ocean data onto a vertical transect. This is useful for
visualizing data along a specific path through the ocean, showing the
vertical structure of ocean properties. These functions allow you to create
cross-sectional plots of temperature, salinity, and other variables.

The following functions are available:

The function :py:func:`mpas_tools.ocean.viz.transect.compute_transect()`
builds a sequence of quads showing the transect intersecting MPAS cells.
This function takes horizontal and vertical mesh information and constructs
a set of quadrilaterals that represent the transect's path through the
MPAS-Ocean mesh. The resulting quads can then be used for plotting or
further analysis.

The remaining functions and those in :ref:`viz_transect_horiz` are lower level
functions that are used by ``compute_transect()`` and are not typically called
directly.

The function
:py:func:`mpas_tools.ocean.viz.transect.find_transect_levels_and_weights()`
constructs a vertical coordinate for a transect and computes interpolation
weights that can be used in ``interp_mpas_to_transect_nodes()`` to performed
linear interpolation.

The function
:py:func:`mpas_tools.ocean.viz.transect.interp_mpas_to_transect_cells()`
interpolates an MPAS-Ocean DataArray to transect cells, keeping constant values
over each MPAS-Ocean cell. This function uses the indices computed by
``find_transect_levels_and_weights()`` to map data from the MPAS-Ocean mesh
onto the transect. The result is an ``xarray.DataArray`` with values
sampled to transect cells.


The function
:py:func:`mpas_tools.ocean.viz.transect.interp_mpas_to_transect_nodes()`
interpolates an MPAS-Ocean DataArray to transect nodes, linearly
interpolating fields between the closest neighboring cells. This function
uses the interpolation weights computed by
``find_transect_levels_and_weights()`` to map data from the MPAS-Ocean mesh
onto the transect. The result is an ``xarray.DataArray`` with values
interpolated to the transect's nodes.

The function
:py:func:`mpas_tools.ocean.viz.transect.interp_transect_grid_to_transect_nodes()`
interpolates a 2D grid of data to transect nodes, linearly interpolating
fields between the closest neighboring cells. This requires that the
``z_transect`` parameter has been passed into function
``find_transect_levels_and_weights()`` (or ``compute_transect()``) uses
weights generated in ``find_transect_levels_and_weights()`` to interpolate data
from the MPAS-Ocean mesh to transect nodes the transect, resulting an
``xarray.DataArray`` with values.
