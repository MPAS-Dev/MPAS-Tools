#!/bin/csh
rm -f locs.dat*
cp centroids.162.dat locs.dat

make clean
make

setenv NAME x1

foreach RES (162 642 2562 10242 40962)
#foreach RES (162 642 2562 10242 40962 163842)

setenv RES 162
echo "&domains" > namelist.input
echo " np = "$RES"" >> namelist.input
echo " locs_as_xyz = .true." >> namelist.input
echo " n_scvt_iterations = 10000" >> namelist.input
cat convergence >> namelist.input
echo "/" >> namelist.input
grid_gen
grid_ref
mv -f grid.nc grid.$NAME.$RES.nc
mv -f locs.dat.out locs.$NAME.$RES.dat
mv -f locs.dat.out.refined locs.dat

end
