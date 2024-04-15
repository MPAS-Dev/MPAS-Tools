cmake_minimum_required (VERSION 3.0.2)
enable_language(Fortran)
project (ocean_cull_wave_mesh)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)

message (STATUS "NETCDF_INCLUDES=${NETCDF_INCLUDES}")
message (STATUS "NETCDF_LIBRARIES=${NETCDF_LIBRARIES}")

include_directories(${NETCDF_INCLUDES})

add_executable (ocean_cull_wave_mesh cull_wave_mesh.F90)
target_link_libraries (ocean_cull_wave_mesh ${NETCDF_LIBRARIES})

install (TARGETS ocean_cull_wave_mesh DESTINATION bin)
