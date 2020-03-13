#!/bin/bash

set -e

source $HOME/miniconda/etc/profile.d/conda.sh
conda activate base
conda build -m "conda_package/travis_ci/linux_${TRAVIS_JOB_NAME}.yaml" "conda_package/recipe"

conda create -y -n test --use-local mpas_tools pytest sphinx mock \
   sphinx_rtd_theme

conda activate test

if [[ "$TRAVIS_BRANCH" == "master" ]]; then
  export DOCS_VERSION="stable"
elif [ -n "$TRAVIS_TAG" ]; then
  # this is a tag build
  export DOCS_VERSION="$TRAVIS_TAG"
else
  DOCS_VERSION=`python -c "import mpas_tools; print(mpas_tools.__version__)"`
  export DOCS_VERSION
fi
cd conda_package/docs || exit 1
make html
