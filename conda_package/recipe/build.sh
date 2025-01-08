#!/usr/bin/env bash

set -x
set -e

cd conda_package
${PYTHON} -m pip install . --no-deps --no-build-isolation -vv

# build and install ocean topography smoothing tool
cd ${SRC_DIR}/ocean/smooth_topo
mkdir build
cd build
cmake \
  -D CMAKE_INSTALL_PREFIX=${PREFIX} \
  -D CMAKE_BUILD_TYPE=Release \
  ..
cmake --build .
cmake --install .

# build and install sea ice partitioning tool
cd ${SRC_DIR}/mesh_tools/seaice_grid_tools
mkdir build
cd build
cmake \
  -D CMAKE_INSTALL_PREFIX=${PREFIX} \
  -D CMAKE_BUILD_TYPE=Release \
  ..
cmake --build .
cmake --install .

# build and install legacy mask creator
cd ${SRC_DIR}/mesh_tools/mesh_conversion_tools
mkdir build
cd build
cmake \
  -D CMAKE_INSTALL_PREFIX=${PREFIX} \
  -D CMAKE_BUILD_TYPE=Release \
  ..
cmake --build .
cp MpasMaskCreator.x ${PREFIX}/bin

# build and install mesh conversion tools
cd ${SRC_DIR}/mesh_tools/mesh_conversion_tools_netcdf_c
mkdir build
cd build
cmake \
  -D CMAKE_INSTALL_PREFIX=${PREFIX} \
  -D CMAKE_BUILD_TYPE=Release \
  ..
cmake --build .
cmake --install .
