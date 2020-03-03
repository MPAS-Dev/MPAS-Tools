#!/bin/bash

check_env () {
  for module in "geometric_features" "mpas_tools" "jigsawpy"
  do
    if python -c "import ${module}"; then
      echo "  ${module} passed"
    else
      echo "  ${module} failed"
      exit 1
    fi
  done

  for exec in "gpmetis" "ffmpeg"
  do
    if ${exec} --help; then
      echo "  ${exec} passed"
    else
      echo "  ${exec} failed"
      exit 1
    fi
  done
}


# Modify the following to choose which e3sm-unified version(s) the python version(s) are installed and whether to make
# an environment with x-windows support under cdat (cdatx) and/or without (nox).  Typically, both environments should
# be created.
versions=(0.1.0)
pythons=(3.7)

default_python=3.7

# Any subsequent commands which fail will cause the shell script to exit
# immediately
set -e

world_read="True"

# The rest of the script should not need to be modified
if [[ $HOSTNAME = "cori"* ]] || [[ $HOSTNAME = "dtn"* ]]; then
  base_path="/global/cfs/cdirs/acme/software/anaconda_envs/base"
  activ_path="/global/cfs/cdirs/acme/software/anaconda_envs"
  group="acme"
elif [[ $HOSTNAME = "acme1"* ]] || [[ $HOSTNAME = "aims4"* ]]; then
  base_path="/usr/local/e3sm_unified/envs/base"
  activ_path="/usr/local/e3sm_unified/envs"
  group="climate"
elif [[ $HOSTNAME = "blueslogin"* ]]; then
  base_path="/lcrc/soft/climate/e3sm-unified/base"
  activ_path="/lcrc/soft/climate/e3sm-unified"
  group="climate"
elif [[ $HOSTNAME = "rhea"* ]]; then
  base_path="/ccs/proj/cli900/sw/rhea/e3sm-unified/base"
  activ_path="/ccs/proj/cli900/sw/rhea/e3sm-unified"
  group="cli900"
elif [[ $HOSTNAME = "cooley"* ]]; then
  base_path="/lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/base"
  activ_path="/lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified"
  group="ccsm"
elif [[ $HOSTNAME = "compy"* ]]; then
  base_path="/compyfs/software/e3sm-unified/base"
  activ_path="/compyfs/software/e3sm-unified"
  group="users"
elif [[ $HOSTNAME = "gr-fe"* ]] || [[ $HOSTNAME = "wf-fe"* ]]; then
  base_path="/usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base"
  activ_path="/usr/projects/climate/SHARED_CLIMATE/anaconda_envs"
  group="climate"
elif [[ $HOSTNAME = "burnham"* ]]; then
  base_path="/home/xylar/Desktop/test_e3sm_unified/base"
  activ_path="/home/xylar/Desktop/test_e3sm_unified"
  group="xylar"
else
  echo "Unknown host name $HOSTNAME.  Add env_path and group for this machine to the script."
  exit 1
fi

if [ ! -d $base_path ]; then
  miniconda=Miniconda3-latest-Linux-x86_64.sh
  wget https://repo.continuum.io/miniconda/$miniconda
  /bin/bash $miniconda -b -p $base_path
  rm $miniconda
fi

# activate the new environment
source ${base_path}/etc/profile.d/conda.sh
conda activate

conda config --add channels conda-forge
conda config --set channel_priority strict
conda update -y --all

for version in "${versions[@]}"
do
  for python in "${pythons[@]}"
  do
    channels="--override-channels -c conda-forge -c defaults -c e3sm"
    packages="python=$python compass=${version}"

    if [[ "$python" == "$default_python" ]]; then
      suffix=""
    else
      suffix="_py${python}"
    fi

    env_name=compass_${version}${suffix}
    if [ ! -d "$base_path/envs/$env_name" ]; then
      echo creating "$env_name"
      conda create -n "$env_name" -y $channels $packages
      conda activate "$env_name"
      conda deactivate
    else
      echo "$env_name" already exists
    fi

    conda activate "$env_name"
    check_env
    conda deactivate

    mkdir -p "$activ_path"

    # make activation scripts
    script=""
    script="${script}"$'\n'"if [ -x \"\$(command -v module)\" ] ; then"
    script="${script}"$'\n'"  module unload python"
    script="${script}"$'\n'"fi"
    script="${script}"$'\n'"source ${base_path}/etc/profile.d/conda.sh"
    script="${script}"$'\n'"conda activate $env_name"
    file_name=$activ_path/load_latest_compass${suffix}.sh
    rm -f "$file_name"
    echo "${script}" > "$file_name"
  done
done

# delete the tarballs and any unused packages
conda clean -y -p -t

# continue if errors happen from here on
set +e

echo "changing permissions on activation scripts"
chown -R "$USER":$group $activ_path/load_latest_compass*
if [ $world_read == "True" ]; then
  chmod -R go+r $activ_path/load_latest_compass*
  chmod -R go-w $activ_path/load_latest_compass*
else
  chmod -R g+r $activ_path/load_latest_compass*
  chmod -R g-w $activ_path/load_latest_compass*
  chmod -R o-rwx $activ_path/load_latest_compass*
fi

echo "changing permissions on environments"
cd $base_path
echo "  changing owner"
chown -R "$USER:$group" .
if [ $world_read == "True" ]; then
  echo "  adding group/world read"
  chmod -R go+rX .
  echo "  removing group/world write"
  chmod -R go-w .
else
  echo "  adding group read"
  chmod -R g+rX .
  echo "  removing group write"
  chmod -R g-w .
  echo "  removing world read/write"
  chmod -R o-rwx .
fi
echo "  done."

