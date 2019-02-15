#!/usr/bin/env bash

set -x
set -e

export NETCDF=$NETCDF_DIR
cd mesh_tools/mesh_conversion_tools
make
