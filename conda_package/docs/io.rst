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

open_dataset and open_mfdataset
===============================

The :py:func:`mpas_tools.io.open_dataset()` and
:py:func:`mpas_tools.io.open_mfdataset()` functions are thin wrappers around
:py:func:`xarray.open_dataset()` and :py:func:`xarray.open_mfdataset()`.  They
select the NetCDF ``engine`` from the module-level
:py:data:`mpas_tools.io.default_engine` variable when an ``engine`` is not
passed explicitly.

This is useful because :py:func:`xarray.open_dataset()` otherwise sniffs the
file for "magic bits" to auto-select a backend, and that probe can crash on
``NETCDF3_64BIT_DATA`` (CDF5) files.  xarray provides no global way to set a
default engine, so ``mpas_tools.io.default_engine`` offers a single,
process-wide knob that applies to both reading and writing.

Example usage:

.. code-block:: python

    import mpas_tools.io
    from mpas_tools.io import open_dataset

    # use the netcdf4 engine everywhere to avoid the CDF5 sniffing crash
    mpas_tools.io.default_engine = 'netcdf4'

    ds = open_dataset('mesh.nc')
