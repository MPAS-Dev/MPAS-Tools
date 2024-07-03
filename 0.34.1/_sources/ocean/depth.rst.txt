.. _ocean_depth:

Adding a Depth Coordinate
=========================

Adding a 1D depth coordinate
----------------------------

The function :py:func:`mpas_tools.ocean.depth.add_depth()` can be used to add
a 1D ``depth`` coordinate that is appropriate for runs with a z-star MPAS-Ocean
mesh.  The coordinate is only approximately correct but is useful for
visualization.

Internally, the depth is computed with
:py:func:`mpas_tools.ocean.depth.compute_depth()`, which could also be called
directly if one has a suitable ``refBottomDepth`` data array indicating the
reference depth of the bottom of each layer in a 1D coordinate that is
independent of both time and horizontal coordinate.


Adding a 3D zMid coordinate
---------------------------

The function :py:func:`mpas_tools.ocean.depth.add_zmid()` can be used to add
a time-independent, 3D ``zMid`` coordinate that is appropriate for runs with
any MPAS-Ocean vertical coordinate that is not a significant function of time.
This is appropriate for both z-star simulations and those with ice-shelf
cavities, which have a more complex vertical coordinate. The coordinate is only
approximately correct because MPAS-Ocean coordinates vary at least slightly
in time (with the undulation of the sea surface).  The time-independent ``zMid``
is appropriate for visualization an analysis that does not need to account for
this time variability.

Internally, the ``zMid`` is computed with
:py:func:`mpas_tools.ocean.depth.compute_zmid()`, which could also be called
directly if one has a suitable ``bottomDepth``, ``maxLevelCell``,
and (reference) ``layerThickness`` data arrays.


Writing a time-dependent, 3D zMid variable
------------------------------------------

The function :py:func:`mpas_tools.ocean.depth.write_time_varying_zmid()` can be
used to write out a time-dependent, 3D ``zMid`` variable to its own file.
This is the "true" MPAS-Ocean vertical coordinate, in contrast to the 1D and
3D time-independent coordinates mentioned above.  However, it requires a
significant amount of disk space so may not be desirable in many contexts.

Internally, the ``<prefix>zMid`` is computed with
:py:func:`mpas_tools.ocean.depth.compute_zmid()` using the time-dependent
``<prefix>layerThickness`` variable, where ``<prefix>`` is a prefix such as
``'timeMonthly_avg_'`` or an empty string (``''``) for no prefix.
