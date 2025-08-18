.. _mesh_interpolation:

.. |---| unicode:: U+2014  .. em dash, trimming surrounding whitespace
   :trim:

*************
Interpolation
*************

Previously, various tools in this package used ``scipy`` for interpolation.
However, the interpolation routines in ``scipy`` are not well suited to
interpolation from regular grids to MPAS meshes—they are slow and very memory
intensive, particularly for large meshes.

For bilinear interpolation from a tensor (regular) lon/lat grid to an MPAS
mesh, it is much faster and more memory-efficient to use the function
:py:func:`mpas_tools.mesh.interpolation.interp_bilin()`.
This function is specifically designed for interpolating from a regular grid
(e.g., longitude/latitude) to the unstructured cell centers of an MPAS mesh,
and should be preferred over generic `scipy` routines for this use case.

Here is an example where we define cell width for an EC mesh (see
:ref:`ec_mesh`), read in longitude and latitude from an MPAS mesh, and
interpolate the cell widths to cell centers on the MPAS mesh.

.. code-block:: python

    import numpy as np
    import netCDF4 as nc4
    from mpas_tools.mesh.interpolation import interp_bilin

    dlon = 1.
    dlat = dlon
    earth_radius = constants['SHR_CONST_REARTH']
    nlon = int(360./dlon) + 1
    nlat = int(180./dlat) + 1
    lon = np.linspace(-180., 180., nlon)
    lat = np.linspace(-90., 90., nlat)

    cellWidth = mdt.EC_CellWidthVsLat(lat)

    # broadcast cellWidth to 2D
    _, cellWidth = np.meshgrid(lon, cellWidthVsLat)

    ds = nc4.Dataset('base_mesh.nc', 'r+')
    lonCell = ds.variables['lonCell'][:]
    latCell = ds.variables['latCell'][:]

    lonCell = np.mod(np.rad2deg(lonCell) + 180., 360.) - 180.
    latCell = np.rad2deg(latCell)

    cellWidthOnMpas = interp_bilin(lon, lat, cellWidth, lonCell, latCell)

.. note::
    - All cell center coordinates must be within the bounds of the input grid.
    - No extrapolation is performed.
    - For geographic data, use degrees for longitude/latitude to avoid
      round-off issues.
