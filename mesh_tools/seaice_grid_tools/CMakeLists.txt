cmake_minimum_required (VERSION 3.5)
enable_language(Fortran)
project (seaice_grid_tools)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)

message (STATUS "NETCDF_INCLUDES=${NETCDF_INCLUDES}")
message (STATUS "NETCDF_LIBRARIES=${NETCDF_LIBRARIES}")

include_directories(${NETCDF_INCLUDES})
add_executable (fix_regrid_output.exe fix_regrid_output.F90)
target_link_libraries (fix_regrid_output.exe ${NETCDF_LIBRARIES})

install (TARGETS fix_regrid_output.exe DESTINATION bin)

