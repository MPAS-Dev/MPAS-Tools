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

versions=(0.1.6)
pythons=(3.8)
mpis=(serial mpich)

default_python=3.8

remove_existing=False

# Any subsequent commands which fail will cause the shell script to exit
# immediately
set -e

world_read="True"
default_mpi=mpich

# The rest of the script should not need to be modified
if [[ $HOSTNAME = "cori"* ]] || [[ $HOSTNAME = "dtn"* ]]; then
  base_path="/global/cfs/cdirs/e3sm/software/anaconda_envs/base"
  activ_path="/global/cfs/cdirs/e3sm/software/anaconda_envs"
  group="e3sm"
  default_mpi=serial
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
  default_mpi=serial
elif [[ $HOSTNAME = "compy"* ]]; then
  base_path="/share/apps/E3SM/conda_envs/base"
  activ_path="/share/apps/E3SM/conda_envs"
  group="users"
  default_mpi=serial
elif [[ $HOSTNAME = "gr-fe"* ]] || [[ $HOSTNAME = "ba-fe"* ]]; then
  base_path="/usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base"
  activ_path="/usr/projects/climate/SHARED_CLIMATE/anaconda_envs"
  group="climate"
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
    for mpi in "${mpis[@]}"
    do
      channels=("--override-channels" "-c" "conda-forge" "-c" "defaults" "-c"
                "e3sm")
      packages=("python=$python")

      if [[ "$python" == "$default_python" ]]; then
        suffix=""
      else
        suffix="_py${python}"
      fi
      if [[ "$mpi" != "$default_mpi" ]]; then
        suffix="${suffix}_${mpi}"
      fi

      if [[ "$mpi" == "mpich" ]]; then
        packages+=("compass=${version}=mpi_mpich*")
      elif [[ "$mpi" == "openmpi" ]]; then
        packages+=("compass=${version}=mpi_openmpi*")
      elif [[ "$mpi" == "serial" ]]; then
        packages+=("compass=${version}=nompi*")
      fi

      env_name=compass_${version}${suffix}
      if [[ "$remove_existing" == "True" ]]; then
        conda remove -y --all -n "$env_name"
      fi

      if [ ! -d "$base_path/envs/$env_name" ]; then
        echo creating "$env_name"
        conda create -n "$env_name" -y "${channels[@]}" "${packages[@]}"
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
      for ext in sh csh
      do
        script=""
        if [[ $ext = "sh" ]]; then
          script="${script}"$'\n'"if [ -x \"\$(command -v module)\" ] ; then"
          script="${script}"$'\n'"  module unload python"
          script="${script}"$'\n'"fi"
        fi
        script="${script}"$'\n'"source ${base_path}/etc/profile.d/conda.${ext}"
        script="${script}"$'\n'"conda activate $env_name"
        file_name=$activ_path/load_latest_compass${suffix}.${ext}
        rm -f "$file_name"
        echo "${script}" > "$file_name"
      done
    done
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

