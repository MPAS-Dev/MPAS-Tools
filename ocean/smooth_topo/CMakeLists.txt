cmake_minimum_required (VERSION 3.5)
enable_language(Fortran)
project (ocean_smooth_topo)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)

message (STATUS "NETCDF_INCLUDES=${NETCDF_INCLUDES}")
message (STATUS "NETCDF_LIBRARIES=${NETCDF_LIBRARIES}")

include_directories(${NETCDF_INCLUDES})

add_executable (ocean_smooth_topo_skip_land_ice smooth_topo_skip_land_ice.F90)
target_link_libraries (ocean_smooth_topo_skip_land_ice ${NETCDF_LIBRARIES})

add_executable (ocean_smooth_topo_before_init smooth_topo_before_init.F90)
target_link_libraries (ocean_smooth_topo_before_init ${NETCDF_LIBRARIES})

add_executable (ocean_smooth_topo_before_init_skip_land_ice smooth_topo_before_init_skip_land_ice.F90)
target_link_libraries (ocean_smooth_topo_before_init_skip_land_ice ${NETCDF_LIBRARIES})

install (TARGETS ocean_smooth_topo_skip_land_ice ocean_smooth_topo_before_init ocean_smooth_topo_before_init_skip_land_ice DESTINATION bin)
