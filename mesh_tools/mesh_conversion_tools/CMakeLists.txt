cmake_minimum_required (VERSION 3.0.2)
project (mesh_conversion_tools)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

include_directories (netcdf-cxx-4.2 .)

set(SOURCES netcdf_utils.cpp netcdf-cxx-4.2/ncvalues.cpp netcdf-cxx-4.2/netcdf.cpp)

add_executable (MpasMeshConverter.x mpas_mesh_converter.cpp ${SOURCES})
target_link_libraries (MpasMeshConverter.x netcdf)

add_executable (MpasCellCuller.x mpas_cell_culler.cpp ${SOURCES})
target_link_libraries (MpasCellCuller.x netcdf)

add_executable (MpasMaskCreator.x mpas_mask_creator.cpp jsoncpp.cpp ${SOURCES})
target_link_libraries (MpasMaskCreator.x netcdf)

install (TARGETS MpasMeshConverter.x MpasCellCuller.x MpasMaskCreator.x DESTINATION bin)
