.. _quick_start:

***********
Quick Start
***********

This guide will help you get started with `mpas_tools` as quickly as possible.

Installing Conda (Miniforge)
============================

If you do not already have conda installed, we recommend installing
`Miniforge <https://github.com/conda-forge/miniforge?tab=readme-ov-file#requirements-and-installers>`_.
This is a lightweight conda installer that default to the conda-forge channel.

To install Miniforge:

1. Download the installer for your platform from the
   `Miniforge releases page <https://github.com/conda-forge/miniforge/releases>`_.

2. Run the installer following the instructions for your operating system.

3. After installation, initialize conda in your shell (if prompted and if
   desired).

Installing mpas_tools from conda-forge
======================================

Once you have conda available, you can install `mpas_tools` into a new
environment:

.. code-block:: bash

    conda create -n mpas_tools_env -c conda-forge mpas_tools
    conda activate mpas_tools_env

This will create and activate a new environment called `mpas_tools_env` with
the latest version of `mpas_tools` and its dependencies.

Verifying the Installation
==========================

To verify that `mpas_tools` is installed correctly, try importing it in Python:

.. code-block:: python

    import mpas_tools
    print(mpas_tools.__version__)

You can also check that command-line tools are available, for example:

.. code-block:: bash

    planar_hex --help

Basic Usage
===========

A simple example of using `mpas_tools` to create a planar hexagonal mesh
that is periodic in one dimension, then cull out the boundary cells (to
remove periodicity in x) and finally to ensure that the mesh has been
converted to the proper MPAS format:

.. code-block:: python

  from mpas_tools.planar_hex import make_planar_hex_mesh
  from mpas_tools.mesh.conversion import convert, cull
  from mpas_tools.io import write_netcdf

  ds = make_planar_hex_mesh(nx=4, ny=4, dc=10e3, nonperiodic_x=True,
                            nonperiodic_y=False)

  ds = cull(ds)
  ds = convert(ds)
  write_netcdf(ds, 'mesh.nc')

The same can be accomplished with the command-line tool:
.. code-block:: bash

    planar_hex --nx 4 --ny 4 --dc 10000 --nonperiodic_x --output base_mesh.nc
    MpasCellCuller.x base_mesh.nc culled_mesh.nc
    MpasMeshConverter.x culled_mesh.nc final_mesh.nc

Further Reading
===============

- :ref:`api_reference`
- :ref:`mesh_creation`
- :ref:`mesh_interpolation`
- :ref:`mesh_interpolation`
