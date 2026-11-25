.. _landice_adjust_bed_to_haf:

************************************
Adjusting Bed to Height Above Flotation
************************************

The ``adjust_bed_to_haf`` command-line tool adjusts bed topography within
grounding line regions to achieve a specified height above flotation (HAF).

Height Above Flotation
======================

Height above flotation (HAF) is a measure of how far an ice column is from
hydrostatic equilibrium with ocean water. It is defined as:

.. math::

    \text{HAF} = z_{\text{surface}} - z_{\text{flotation}}

where :math:`z_{\text{surface}}` is the ice surface elevation and
:math:`z_{\text{flotation}}` is the surface elevation at which the ice would
be in flotation equilibrium.

For ice in flotation equilibrium:

.. math::

    \rho_{\text{ice}} \cdot H = \rho_{\text{ocean}} \cdot d

where :math:`H` is ice thickness, :math:`d` is the draft (depth below sea
level), and :math:`\rho_{\text{ice}}` and :math:`\rho_{\text{ocean}}` are ice
and ocean water densities, respectively.

Usage
=====

Basic usage:

.. code-block:: bash

    adjust_bed_to_haf \
        --mesh input_mesh.nc \
        --geojson grounding_line.geojson \
        --output output_mesh.nc \
        --target-haf 10.0

This will adjust the bed topography in cells within the grounding line polygons
to achieve a height above flotation of 10 meters.

Command-Line Options
====================

Required Arguments
------------------

``-m, --mesh MESH_FILE``
    MALI mesh file in NetCDF format containing ice thickness and bed topography

``-g, --geojson GEOJSON_FILE``
    GeoJSON file containing grounding line delineations as polygon geometries

Optional Arguments
------------------

``-o, --output OUTPUT_FILE``
    Output mesh file. If not specified, the input mesh file is modified in place.

``--target-haf HAF``
    Target height above flotation in meters (default: 10.0)

``--sea-level LEVEL``
    Sea level in meters (default: 0.0)

``--rho-ice DENSITY``
    Ice density in kg/m³ (default: 910.0)

``--rho-ocean DENSITY``
    Ocean water density in kg/m³ (default: 1028.0)

``--thickness-var VARNAME``
    Name of thickness variable in mesh file (default: 'thickness')

``--bed-var VARNAME``
    Name of bed topography variable in mesh file (default: 'bedTopography')

Example
=======

The following example adjusts bed topography within Thwaites Glacier grounding
line regions to achieve a 15-meter height above flotation:

.. code-block:: bash

    adjust_bed_to_haf \
        --mesh AIS_4to20km_mesh.nc \
        --geojson Thwaites_GL_2014_pinning_points.geojson \
        --output AIS_4to20km_mesh_adjusted.nc \
        --target-haf 15.0 \
        --rho-ice 918.0 \
        --rho-ocean 1028.0

The tool will:

1. Load the grounding line polygons from the GeoJSON file
2. Identify all mesh cells that fall within these polygons
3. Calculate the current height above flotation in these regions
4. Adjust the bed topography to achieve the target HAF while preserving ice thickness
5. Update the output file with metadata about the changes

Technical Details
=================

The tool calculates the required bed elevation :math:`b` to achieve a target
HAF given ice thickness :math:`H`:

.. math::

    b = z_{\text{sea}} + \text{HAF}_{\text{target}} - H \cdot \frac{\rho_{\text{ice}}}{\rho_{\text{ocean}}}

This ensures that with the existing ice thickness, the ice surface will be at
the specified height above the flotation level.

Notes
=====

- The GeoJSON file should use WGS 84 (EPSG:4326) coordinates (longitude, latitude)
- The tool converts mesh coordinates from radians to degrees for comparison with GeoJSON
- If no cells are found within the grounding line polygons, check that the coordinate systems are compatible
- The tool updates the file's ``history`` and ``comment`` global attributes to document the modifications
- Ice thickness is not modified; only bed topography is adjusted

Dependencies
============

Required Python packages:

- netCDF4
- numpy
- shapely

Optional (for better performance):

- geopandas (falls back to json + shapely if not available)

References
==========

Wild, C. T., Alley, K. E., Muto, A., Pettit, E. C., Scambos, T. A., & Truffer, M. (2022).
Thwaites Glacier 2014 and 2019/20 grounding line positions.
U.S. Antarctic Program (USAP) Data Center.
doi: 10.15784/601499
