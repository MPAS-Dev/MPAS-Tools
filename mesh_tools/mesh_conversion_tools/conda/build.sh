#!/usr/bin/env bash

set -x
set -e

cd mesh_tools/mesh_conversion_tools

export CC=${GXX}
export CFLAGS="-O3 -std=c++0x -fopenmp -lstdc++"

make

install -d ${PREFIX}/bin/
for exec in MpasMeshConverter.x MpasCellCuller.x MpasMaskCreator.x mark_horns_for_culling.py
do
  install -m 755 ${exec} ${PREFIX}/bin/
done