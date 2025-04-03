cmake_minimum_required (VERSION 3.5)
project (mesh_conversion_tools)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable (MpasMeshConverter.x mpas_mesh_converter.cpp)
target_link_libraries (MpasMeshConverter.x netcdf)

add_executable (MpasCellCuller.x mpas_cell_culler.cpp)
target_link_libraries (MpasCellCuller.x netcdf)

#add_executable (MpasMaskCreator.x mpas_mask_creator.cpp jsoncpp.cpp)
#target_link_libraries (MpasMaskCreator.x netcdf)

install (TARGETS MpasMeshConverter.x MpasCellCuller.x DESTINATION bin)
