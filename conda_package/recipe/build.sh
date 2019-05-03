#!/usr/bin/env bash

set -x
set -e

cp -r ocean landice visualization mesh_tools conda_package

cd conda_package
${PYTHON} -m pip install . --no-deps -vv

cd mesh_tools/mesh_conversion_tools

export CXX=${GXX}
export CFLAGS="-O3 -std=c++0x -fopenmp -lstdc++"

make

install -d ${PREFIX}/bin/
for exec in MpasMeshConverter.x MpasCellCuller.x MpasMaskCreator.x
do
  install -m 755 ${exec} ${PREFIX}/bin/
done

