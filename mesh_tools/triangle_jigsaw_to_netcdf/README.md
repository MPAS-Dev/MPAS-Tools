# Mesh Creation and Conversion Steps

## Mesh creation examples
Example scripts in the `examples` folder can be run via `./build_mesh.sh
examples/mpas_uniform.m` for example, with output in
`examples/mpas_uniform.nc`, for example.  Note, paths need to be updated
at top of `build_mesh.sh` as well as in `examples/jigsaw_path_locations.m`.
Currently, `jigsaw-geo-matlab` is assumed to be in the same directory as `MPAS-Tools`.

## General JIGSAW mesh conversion instructions
1. Produce a JIGSAW mesh, e.g., example.msh, from https://github.com/dengwirda/jigsaw-geo-matlab
2. `./triangle_jigsaw_to_netcdf.py -m example.msh -s`
3. `./MpasMeshConverter.x grid.nc`
4. Final mesh `mesh.nc` can then be used to create our initial condition files.

## General TRIANGLE mesh conversion instructions
1. Produce a TRIANGLE mesh, e.g., produced from http://www.netlib.org/voronoi/triangle.zip
2. `./triangle -p example.poly`
3. `./triangle_jigsaw_to_netcdf.py -n example.node -e example.ele`
4. `./MpasMeshConverter.x grid.nc`
5. Final mesh `mesh.nc` can then be used to create our initial condition files.
