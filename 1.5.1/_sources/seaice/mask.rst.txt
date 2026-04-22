.. _seaice_mask:

Mask
====

The :py:mod:`mpas_tools.seaice.mask` module contains a function for
manipulating sea-ice region masks.

.. _seaice_mask_extend_seaice_mask:

Extending a Mask
----------------

The function :py:func:`mpas_tools.seaice.mask.extend_seaice_mask()` is used to
extend a sea-ice "presence" mask that covers the area where sea-ice is present
by a given distance.  This is useful as part of creating a sea-ice
:ref:`seaice_partition`.
