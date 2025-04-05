import argparse
import os
import platform
import subprocess
import time

import numpy

from mpas_tools.logging import check_call


def jigsaw_driver(
    cellWidth,
    x,
    y,
    on_sphere=True,
    earth_radius=6371.0e3,
    geom_points=None,
    geom_edges=None,
    logger=None,
):
    """
    A function for building a jigsaw mesh

    Parameters
    ----------
    cellWidth : ndarray
        The size of each cell in the resulting mesh as a function of space

    x, y : ndarray
        The x and y coordinates of each point in the cellWidth array (lon and
        lat for spherical mesh)

    on_sphere : logical, optional
        Whether this mesh is spherical or planar

    earth_radius : float, optional
        Earth radius in meters

    geom_points : ndarray, optional
        list of point coordinates for bounding polygon for planar mesh

    geom_edges : ndarray, optional
        list of edges between points in geom_points that define the bounding
        polygon

    logger : logging.Logger, optional
        A logger for the output if not stdout
    """
    try:
        import jigsawpy
        from jigsawpy.savejig import savejig
    except ImportError as err:
        raise ImportError(
            'JIGSAW and/or jigsaw-python is not installed. '
            'Please install them in the conda environment.'
        ) from err
    # Authors
    # -------
    # Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    # setup files for JIGSAW
    opts = jigsawpy.jigsaw_jig_t()
    opts.geom_file = 'mesh.msh'
    opts.jcfg_file = 'mesh.jig'
    opts.mesh_file = 'mesh-MESH.msh'
    opts.hfun_file = 'mesh-HFUN.msh'

    # save HFUN data to file
    hmat = jigsawpy.jigsaw_msh_t()
    if on_sphere:
        hmat.mshID = 'ELLIPSOID-GRID'
        hmat.xgrid = numpy.radians(x)
        hmat.ygrid = numpy.radians(y)
    else:
        hmat.mshID = 'EUCLIDEAN-GRID'
        hmat.xgrid = x
        hmat.ygrid = y
    hmat.value = cellWidth
    jigsawpy.savemsh(opts.hfun_file, hmat)

    # define JIGSAW geometry
    geom = jigsawpy.jigsaw_msh_t()
    if on_sphere:
        geom.mshID = 'ELLIPSOID-MESH'
        geom.radii = earth_radius * 1e-3 * numpy.ones(3, float)
    else:
        geom.mshID = 'EUCLIDEAN-MESH'
        geom.vert2 = geom_points
        geom.edge2 = geom_edges
    jigsawpy.savemsh(opts.geom_file, geom)

    # build mesh via JIGSAW!
    opts.hfun_scal = 'absolute'
    opts.hfun_hmax = float('inf')
    opts.hfun_hmin = 0.0
    opts.mesh_dims = +2  # 2-dim. simplexes
    opts.optm_qlim = 0.9375
    opts.verbosity = +1

    savejig(opts.jcfg_file, opts)
    check_call(['jigsaw', opts.jcfg_file], logger=logger)


def build_jigsaw(logger=None, clone=False):
    """
    Build the JIGSAW and JIGSAW-Python tools using conda-forge compilers

    Parameters
    ----------
    logger : logging.Logger, optional
        A logger for the output if not stdout

    clone : bool, optional
        If True, clone the jigsaw-python repository from github
        and build it. If False, just build the existing repository.
    """
    conda_env_path = os.getenv('CONDA_PREFIX')
    if conda_env_path is None:
        raise EnvironmentError(
            'The CONDA_PREFIX environment variable is not defined. '
            'Please activate a conda environment where you want to install '
            'jigsaw and jigsawpy before running this function.'
        )

    conda_exe = os.getenv('CONDA_EXE')
    if conda_exe is None:
        raise EnvironmentError(
            'The CONDA_EXE environment variable is not defined. '
            'Please ensure conda is installed and accessible.'
        )
    conda_base = os.path.dirname(os.path.dirname(conda_exe))

    conda_sh_path = os.path.join(conda_base, 'etc', 'profile.d', 'conda.sh')
    conda_env_name = os.getenv('CONDA_DEFAULT_ENV')
    if conda_env_name is None:
        raise EnvironmentError(
            'The CONDA_DEFAULT_ENV environment variable is not defined. '
            'Please ensure a conda environment is activated.'
        )

    activate_env = (
        f'source {conda_sh_path} && conda activate {conda_env_name} && '
    )

    # remove conda jigsaw and jigsaw-python
    t0 = time.time()
    commands = f'{activate_env}conda remove -y --force-remove jigsaw jigsawpy'
    try:
        check_call(commands, logger=logger, executable='/bin/bash', shell=True)
    except subprocess.CalledProcessError:
        # ignore errors if not installed
        pass

    if clone:
        commands = (
            f'{activate_env}'
            f'rm -rf jigsaw-python && '
            f'git clone https://github.com/dengwirda/jigsaw-python.git'
        )
        check_call(commands, logger=logger, executable='/bin/bash', shell=True)

    # add build tools to deployment env, not polaris env
    jigsaw_build_deps = (
        'cxx-compiler cmake make libnetcdf setuptools numpy scipy'
    )
    if platform.system() == 'Linux':
        jigsaw_build_deps = f'{jigsaw_build_deps} sysroot_linux-64=2.17'
        netcdf_lib = f'{conda_env_path}/lib/libnetcdf.so'
    elif platform.system() == 'Darwin':
        jigsaw_build_deps = (
            f'{jigsaw_build_deps} macosx_deployment_target_osx-64=10.13'
        )
        netcdf_lib = f'{conda_env_path}/lib/libnetcdf.dylib'

    cmake_args = f'-DCMAKE_BUILD_TYPE=Release -DNETCDF_LIBRARY={netcdf_lib}'

    print('Install dependencies\n')
    # Install dependencies
    commands = f'{activate_env}conda install -y {jigsaw_build_deps}'
    check_call(commands, logger=logger, executable='/bin/bash', shell=True)

    print('Building JIGSAW\n')
    # Build JIGSAW
    commands = (
        f'{activate_env}'
        f'cd jigsaw-python/external/jigsaw && '
        f'rm -rf tmp && '
        f'mkdir tmp && '
        f'cd tmp && '
        f'cmake .. {cmake_args} && '
        f'cmake --build . --config Release --target install --parallel 4'
    )
    check_call(commands, logger=logger, executable='/bin/bash', shell=True)

    print('Installing JIGSAW into JIGSAW-Python\n')
    # Set up JIGSAW-Python
    commands = (
        f'{activate_env}'
        f'cd jigsaw-python && '
        f'rm -rf jigsawpy/_bin jigsawpy/_lib && '
        f'cp -r external/jigsaw/bin/ jigsawpy/_bin && '
        f'cp -r external/jigsaw/lib/ jigsawpy/_lib'
    )
    check_call(commands, logger=logger, executable='/bin/bash', shell=True)

    print('Installing JIGSAW-Python\n')
    commands = (
        f'{activate_env}'
        f'cd jigsaw-python && '
        f'python -m pip install --no-deps --no-build-isolation -e . && '
        f'cp jigsawpy/_bin/* ${{CONDA_PREFIX}}/bin'
    )
    check_call(commands, logger=logger, executable='/bin/bash', shell=True)

    t1 = time.time()
    total = int(t1 - t0 + 0.5)
    message = f'JIGSAW install took {total:.1f} s.'
    if logger is None:
        print(message)
    else:
        logger.info(message)


def main_build_jigsaw():
    """
    Entry point for building JIGSAW and JIGSAW-Python tools.
    """

    parser = argparse.ArgumentParser(
        description='Build JIGSAW and JIGSAW-Python tools.'
    )
    parser.add_argument(
        '--clone',
        dest='clone',
        action='store_true',
        help=(
            'Clone the jigsaw-python repository into the jigsaw-python dir '
            'for you.  Otherwise, you must already have cloned it.'
        ),
    )
    args = parser.parse_args()

    build_jigsaw(clone=args.clone)
