.. _ocean_streamfunction:

Computing streamfunctions
=========================

Computing the barotropic streamfunction
---------------------------------------

The function :py:func:`mpas_tools.ocean.compute_barotropic_streamfunction()`
computes the barotproic streamfunction at vertices on the MPAS-Ocean grid.
The function takes a dataset containing an MPAS-Ocean mesh and another with
``normalVelocity`` and ``layerThickness`` variables (possibly with a
``timeMonthly_avg_`` prefix).  The streamfunction is computed only over the
range of (positive-down) depths provided and at the given time index.
