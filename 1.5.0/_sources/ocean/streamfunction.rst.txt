.. _ocean_streamfunction:

Computing streamfunctions
=========================

Computing the barotropic streamfunction
---------------------------------------

The function :py:func:`mpas_tools.ocean.compute_barotropic_streamfunction()`
computes the barotropic streamfunction at vertices on the MPAS-Ocean grid.
The function takes a dataset containing an MPAS-Ocean mesh and another with
``normalVelocity`` and ``layerThickness`` variables (possibly with a
``timeMonthly_avg_`` prefix). The streamfunction is computed only over the
range of (positive-down) depths provided and at the given time index.

Optionally, the Gent-McWilliams bolus velocity (``normalGMBolusVelocity``) and
the submesoscale velocity (``normalMLEvelocity``) can be included in the
vertically integrated velocity calculation by setting the corresponding
arguments ``include_bolus`` and ``include_submesoscale`` to ``True``.

For large meshes, performance and memory usage can be improved by specifying
the ``horiz_chunk`` argument to control the number of edges processed at once,
and by providing a ``tmp_dir`` argument to use a temporary directory for
intermediate files.

Shifting the barotropic streamfunction
--------------------------------------

The function :py:func:`mpas_tools.ocean.shift_barotropic_streamfunction()`
can be used to shift the barotropic streamfunction so that its mean value is
zero on the boundary within a specified latitude range. This is useful for
setting a consistent reference for the streamfunction, especially when
comparing results across different regions or simulations.

The function takes as input the barotropic streamfunction on vertices,
a latitude range (in degrees), the ``cellsOnVertex`` and ``latVertex``
arrays from the mesh, and optionally a logger. It returns a shifted
streamfunction with the mean value on the boundary (within the latitude
range) subtracted.
