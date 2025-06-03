.. _io:

*********
I/O Tools
*********

The :py:mod:`mpas_tools.io` module provides utilities for reading and writing
NetCDF files, especially for compatibility with MPAS mesh and data conventions.

write_netcdf
============

The :py:func:`mpas_tools.io.write_netcdf()` function writes an
``xarray.Dataset`` to a NetCDF file, ensuring MPAS compatibility (e.g.,
converting int64 to int32, handling fill values, and updating the history
attribute). It also supports writing in various NetCDF formats, including
conversion to ``NETCDF3_64BIT_DATA`` using ``ncks`` if needed.

Example usage:

.. code-block:: python

    import xarray as xr
    from mpas_tools.io import write_netcdf

    # Create a simple dataset
    ds = xr.Dataset({'foo': (('x',), [1, 2, 3])})
    write_netcdf(ds, 'output.nc')
