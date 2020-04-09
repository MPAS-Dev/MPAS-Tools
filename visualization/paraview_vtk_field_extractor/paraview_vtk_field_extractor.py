#!/usr/bin/env python
"""
Name: paraview_vtk_field_extractor.py
Authors: Doug Jacobsen, Xylar Asay-Davis

This script is used to extract a field from a time series of NetCDF files as
VTK files for plotting in paraview.

It can extract a field across multiple files by passing in a regular expression
for the filename patter. As an example, one can run the script using:

    ./paraview_vtk_field_extractor.py -v areaCell,latVertex -f "hist.comp.*.nc"

To extract a time series of areaCell,latVertex that spans multiple files.
By default, time-independent fields on cells are written to a file
vtk_files/staticFieldsOnCells.vtp
and time-dependent fields on cells are written to
vtk_files/timeDependentFieldsOnCells.pvd
vtk_files/time_series/timeDependentFieldsOnCells.0.vtp
vtk_files/time_series/timeDependentFieldsOnCells.1.vtp
...
and similarly for edges and vertices.  Time-independent fields can be
included in each time step of the time-dependent fields for with
the --combine flag.  This allows time-dependent and -independent fields
to be combined in filters within Paraview at the expense of considerable
additional storage space.

Indices for extra dimensions can either be supplied at runtime at a prompt
or through the -d flag.  For each extra dimension, the use can specifiy
nothing (an empty string, meaning skip any fields with this dimension),
a single index, or a comma-separated list of indices or a range of indices
indices (separated by 1 or 2 colons).  For example,

    -d maxEdges= nVertLeves=0:10:2 nParticles=0,2,4,6,8

will ignore any fields with dimension maxEdges, extract every other layer from
the first 10 vertical levels (each into its own field) and extract the five
specified particles.

An index array can also be specified in this way (and these can be mixed with
integer indices in a comma-separated list but not in a colon-separated range):

    -d nVertLeves=0,maxLevelCell

will extract fields from the first vertical level and the vertical level with
index given by maxLevelCell.

The extractor includes optional support for extracting geometry appropriate
for displaying variables at the depth of a topographic feature (typically the
top or bottom of the domain) for MPAS components with a spatially variable
top or bottom index (e.g. `maxLevelCell` in MPAS-Ocean).  This is accomplished
with flags such as:

    --topo_dim=nVertLevels --topo_cell_index=maxLevelCell

Fields on cells are sampled at the topographic index and the geometry includes
polygons corresponding to edges so that vertical faces between adjacent cells
can be displayed.  Fields are extracted as normal except that they are sampled
as point data rather than cell data, allowing computations in ParaView to
display the topography.  A mask field is also included indicating which parts
of edge polygons correspond to the boundary of the domain (boundaryMask == 1)
and which parts of cell and edge polygons are interior (boundaryMask == 0).
Together, this can be used to plot topography by using a calculator filter like
the following:

    coords*(1.0 + 100.0/mag(coords)*((1 - boundaryMask)*(-bottomDepth)
                                     + 10.0*boundaryMask))

If this is entered into a Calculator Filter in ParaView with the "coordinate
result" box checked, the result will to display the MPAS-Ocean topography,
exaggerated by a factor of 100, with a value equivalent to 10 m along boundary
points of edge polygons (a "water-tight" surface).
"""

from mpas_tools.viz.paraview_extractor import main


if __name__ == "__main__":
    main()

