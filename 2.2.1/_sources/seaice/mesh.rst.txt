.. _seaice_mesh:

Mesh
====

The :py:mod:`mpas_tools.seaice.mesh` module contains several functions for
creating scrip files for sea-ice meshes and regular grids.

.. _seaice_mesh_mpas_scrip:

MPAS-Seaice SCRIP files
-----------------------

The functions :py:func:`mpas_tools.seaice.mesh.make_mpas_scripfile_on_cells()`
and :py:func:`mpas_tools.seaice.mesh.make_mpas_scripfile_on_vertices()` are for
creating scrip files for fields on cells and vertices, respectively.

These are both created using a lower-level function
:py:func:`mpas_tools.seaice.mesh.write_scrip_file()`.


.. _seaice_mesh_2d_scrip:

SCRIP files for 2D grids
------------------------

The function :py:func:`mpas_tools.seaice.mesh.write_2D_scripfile()` is for
creating scrip files for a regular, 2D latitude/longitude grid.
