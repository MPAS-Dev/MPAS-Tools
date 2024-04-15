cmake_minimum_required (VERSION 3.0.2)
enable_language(Fortran)
project (ocean_cull_wave_mesh)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)

message (STATUS "NETCDF_INCLUDES=${NETCDF_INCLUDES}")
message (STATUS "NETCDF_LIBRARIES=${NETCDF_LIBRARIES}")

include_directories(${NETCDF_INCLUDES})

add_executable (ocean_cull_wave_mesh cull_wave_mesh.f90 in_cell.f90 edge_connectivity_mod.f90 fix_elements.f90 globals.f90 grid_file_mod.f90 kdtree2.f90 read_write_gmsh.f90 write_vtk.f90)
target_link_libraries (ocean_cull_wave_mesh ${NETCDF_LIBRARIES})

install (TARGETS ocean_cull_wave_mesh DESTINATION bin)
