#!/usr/bin/env bash
# Directs process to build MPAS mesh using JIGSAW
# Phillip J Wolfram, 01/19/2018

MATLAB=/Applications/MATLAB_R2015b.app/bin/matlab
MESHCONVERTER=../mesh_conversion_tools/MpasMeshConverter.x
JIGSAW2NETCDF=./
NAME=${1%.*}

echo 'Starting workflow to build grid from '$NAME.m':'

echo 'Build mesh using JIGSAW ...'
CMD='try, run('"'$NAME.m'"'), catch, exit(1), end, exit(0);'
$MATLAB -nodesktop -nodisplay -nosplash  -r  "$CMD"
echo 'done'

echo 'Convert to netcdf file (grid.nc) ...'
python $JIGSAW2NETCDF/triangle_jigsaw_to_netcdf.py -m ${NAME}-MESH.msh -s
echo 'done'

echo 'Convert to MPAS mesh (mesh.nc) ...'
$MESHCONVERTER grid.nc
echo 'done'

echo 'Removing grid.nc and renaming mesh.nc to '$NAME'.nc ...'
rm grid.nc
mv mesh.nc $NAME.nc
echo 'done'
