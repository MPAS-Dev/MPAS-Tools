.. _seaice_regions:

Region masks
============

The :py:mod:`mpas_tools.seaice.regions` module contains a function for
creating masks to help with graph partitioning.

.. _seaice_regions_make_regions_file:

Make a region mask for partitioning
-----------------------------------

The function :py:func:`mpas_tools.seaice.regions.make_regions_file()` is used
to create a ``region`` field with different integer values for different
regions that are used as part of creating a sea-ice :ref:`seaice_partition`.
