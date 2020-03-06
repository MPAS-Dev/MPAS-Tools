#!/bin/bash

conda build -c conda-forge -c e3sm .

#anaconda upload -u e3sm ${HOME}/miniconda3/conda-bld/noarch/compass*.tar.bz2
