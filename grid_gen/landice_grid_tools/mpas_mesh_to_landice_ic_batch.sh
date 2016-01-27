#!/bin/bash
#MSUB -l walltime=02:00:00
#MSUB -l nodes=1
#MSUB -j oe

source ~/setup_gnu_env.sh

# This script is a template for a batch script that will perform the various steps required
# to go from a base MPAS mesh to a culled landice mesh with initial conditions.
# It's main purpose is to be used on clusters with a batch system for situations
# where the larger memory available on those machines is required.
# In addition to setting the variables below, it is recommended to look over all the steps
# and ensure they are happening in the proper order and with the correct flags for
# your particular situation.


# Location of MPAS-Tools repo
TOOLS=/users/mhoffman/mpas/MPAS-Tools

# Command used on the machine for running an executable
RUNCMD="mpirun -n 1 "

# The name of the mesh with which processing
STARTMESH=mpas.nc

# The name of the initial condition file from which initial conditions will be interpolated
INTERPFILE=/usr/projects/climate/mhoffman/AIS_IC_data/Antarctica-1km.BISICLES.CISM-style.nc

# Amount to translate the mesh in x and y
SHIFTX=0.0
SHIFTY=0.0

# Projection to use for adding lat/lon values.  See that script for options or to add new ones
PROJ=ais-bedmap2

# ==========================

date

# shift by some amount
#$RUNCMD $TOOLS/grid_gen/planar_grid_transformations/translate_planar_grid.py -f $STARTMESH -x $SHIFTX -y $SHIFTY 
date

# generate culling mask
$RUNCMD $TOOLS/grid_gen/landice_grid_tools/define_cullMask.py -f $STARTMESH -m 4  # method 4 removed junk from near the edges
date

# Cull mesh
$RUNCMD $TOOLS/grid_gen/mesh_conversion_tools/MpasCellCuller.x $STARTMESH culled1.nc
date

# Add lat/lon
$RUNCMD $TOOLS/grid_gen/planar_grid_transformations/set_lat_lon_fields_in_planar_grid.py -f culled1.nc -p $PROJ
date

# convert to LI mesh - no reason to have many vertical levels yet since we are going to cull
$RUNCMD $TOOLS/grid_gen/landice_grid_tools/create_landice_grid_from_generic_MPAS_grid.py -i culled1.nc -o landice_grid_full.nc -l 2 --beta
date

# Interpolate thickness
$RUNCMD $TOOLS/grid_gen/landice_grid_tools/interpolate_cism_grid_to_mpas_grid.py -m landice_grid_full.nc -c $INTERPFILE --thickness-only
date

# generate culling mask
$RUNCMD $TOOLS/grid_gen/landice_grid_tools/define_cullMask.py -f landice_grid_full.nc -m 2
date

# Cull mesh
$RUNCMD $TOOLS/grid_gen/mesh_conversion_tools/MpasCellCuller.x landice_grid_full.nc culled2.nc
date

# Make LI mesh
$RUNCMD $TOOLS/grid_gen/landice_grid_tools/create_landice_grid_from_generic_MPAS_grid.py -i culled2.nc -o landice_grid.nc -l 10 --beta
date

# Interpolate everything
$RUNCMD $TOOLS/grid_gen/landice_grid_tools/interpolate_cism_grid_to_mpas_grid.py -m landice_grid.nc -c $INTERPFILE 
date

echo "All done."
