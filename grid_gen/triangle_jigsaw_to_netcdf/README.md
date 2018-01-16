# Mesh Conversion Steps

## JIGSAW mesh
1. Produce a JIGSAW mesh, e.g., example.msh, from https://github.com/dengwirda/jigsaw-geo-matlab
2. `./triangle_jigsaw_to_netcdf.py -m example.msh -s`
3. `./MpasMeshConverter.x grid.nc`
4. Final mesh mesh.nc can then be used to create our initial condition files.

## TRIANGLE mesh
1. Produce a TRIANGLE mesh, e.g., produced from http://www.netlib.org/voronoi/triangle.zip
2. `./triangle -p example.poly`
3. `./triangle_jigsaw_to_netcdf.py -n example.node -e example.ele`
4. `./MpasMeshConverter.x grid.nc`
5. Final mesh mesh.nc can then be used to create our initial condition files.
