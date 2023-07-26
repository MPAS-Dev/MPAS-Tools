.. _seaice_partitions:

Graph partitioning
==================

The :py:mod:`mpas_tools.seaice.partition` module contains a function for
creating graph partition files for MPAS-Seaice that are better balanced than
those created from Metis tools directly.

.. _seaice_partitions_:

Running from compass
--------------------

One way to run the tools is from compass using the
`files_for_e3sm <https://mpas-dev.github.io/compass/latest/users_guide/ocean/test_groups/global_ocean.html#files-for-e3sm-for-an-existing-mesh>`_
test case.

This has the advantage that it can run with a version of ESMF that has been
compiled with system compilers for compass.  Compass also automatically
downloads and links the files needed for determining the regions of sea-ice
coverage.  However, if you are unfamiliar with compass, there may be a learning
curve involved in setting up and running the test case.

Conda environment
-----------------

The other preferred way to use the sea ice partitioning tool is through the
mpas_tools conda package.  To install it, first install
`Mambaforge <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh>`_
(if you don't already have Miniconda3):

To activate it, run:

.. code-block:: bash

    source ~/mambaforge/etc/profile.d/conda.sh
    source ~/mambaforge/etc/profile.d/mamba.sh

To create a new conda environment for ``mpas_tools``, run:

.. code-block:: bash

    mamba activate
    mamba create -y -n mpas_tools python=3.11 mpas_tools "esmf=*=nompi*"

This will create a new conda environment called ``mpas_tools`` that contains
the ``mpas_tools`` package and also the version of ESMF without MPI support.
This is necessary because the version with MPI support (the default) doesn't
typically work on HPC.

Each time you want to run the sea-ice partitioning tools, run:

.. code-block:: bash

    source ~/mambaforge/etc/profile.d/conda.sh
    source ~/mambaforge/etc/profile.d/mamba.sh
    mamba activate mpas_tools

All the tools (including ``fix_regrid_output.exe`` built from the Fortran code)
are part of the ``mpas_tools`` conda package.

You will also need an MPAS mesh file to partition.  You do not need to pass
a location for the MPAS cell culler (``-c``) or Metis (``-g``) because theses
will be found automatically in the conda package


Graph partition tools
---------------------
The tools ``prepare_seaice_partitions`` and ``create_seaice_partitions`` are
used to create MPAS-Seaice graph partitions that are better balanced so that
each processor "owns" cells from both polar and equatorial regions.

.. code-block:: none

    $ prepare_seaice_partitions --help
    usage: prepare_seaice_partitions [-h] -i MESHFILENAMESRC -p FILENAMEDATA -m
                                     MESHFILENAMEDST -o OUTPUTDIR

    Perform preparatory work for making seaice partitions.

    options:
      -h, --help            show this help message and exit
      -i MESHFILENAMESRC, --inputmesh MESHFILENAMESRC
                            MPAS mesh file for source regridding mesh.
      -p FILENAMEDATA, --presence FILENAMEDATA
                            Input ice presence file for source mesh.
      -m MESHFILENAMEDST, --outputmesh MESHFILENAMEDST
                            MPAS mesh file for destination regridding mesh.
      -o OUTPUTDIR, --outputDir OUTPUTDIR
                            Output direct

The input mesh file contains the MPAS-Seaice mesh fields for the mesh used to
create the "presence" file.  The "presence" file itself contains an
``icePresence`` field that indicates where sea ice might be present. We
typically use a
`60-km mesh file <https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/partition/seaice_QU60km_polar.nc>`_
and the corresponding
`presence file <https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/partition/icePresent_QU60km_polar.nc>`_.
The presence file will be regridded to the given output MPAS-Seaice mesh. Here the ice
presence was determined as any cell in the mesh to have had ice present at any time
during a 50 year MPAS-Seaice standalone simulation with the above 60-km mesh file.
The output directory is often the current directory.

After this preprocessing has finished, the ``create_seaice_partitions`` tool
can be run one or more times.  It is significantly more efficient to provide
a list of processor numbers than to call the tool for each processor number
separately.

.. code-block:: none

    $ create_seaice_partitions --help
    usage: create_seaice_partitions [-h] -m MESHFILENAME -o OUTPUTDIR
                                    [-c MPASCULLERLOCATION] [-p OUTPUTPREFIX] [-x]
                                    [-g METIS] [-n NPROCS] [-f NPROCSFILE]

    Create sea-ice partitions.

    options:
      -h, --help            show this help message and exit
      -m MESHFILENAME, --outputmesh MESHFILENAME
                            MPAS mesh file for destination regridding mesh.
      -o OUTPUTDIR, --outputDir OUTPUTDIR
                            Output directory for temporary files and partition
                            files.
      -c MPASCULLERLOCATION, --cullerDir MPASCULLERLOCATION
                            Location of MPAS MpasCellCuller.x executable.
      -p OUTPUTPREFIX, --prefix OUTPUTPREFIX
                            prefix for output partition filenames.
      -x, --plotting        create diagnostic plotting file of partitions
      -g METIS, --metis METIS
                            name of metis utility
      -n NPROCS, --nProcs NPROCS
                            number of processors to create partition for.
      -f NPROCSFILE, --nProcsFile NPROCSFILE
                            number of processors to create partition for.

The mesh filename provides the desired MPAS-Seaice mesh, the same as the
destination mesh for ``prepare_seaice_partitions``.  The output directory
is often the current directory.  A directory containing the
``MpasCellCuller.x`` tool can be provided but by default it will be found in
your path as part of the ``mpas_tools`` conda package.  The output prefix will
be prepended onto each graph partition file, and defaults to ``graph.info``.
The Metis tool is nearly always ``gpmetis``, the default, and must be available
in your path (which is the case if you use ``mpas_tools`` conda package).
One graph partition file is created for each number of processors (one or more
integers) provided.  Alternatively, these can be listed, one value on each
line, in a file. You can optionally save a NetCDF file with partition
information ``partition_diag.nc``, which will contain a ``partition_{nProcs}``
field for each number of processors requested.

A simplified tool, primarily intended for use on LCRC machines Anvil and
Chrysalis, has only a few arguments:

.. code-block:: none

    $ simple_seaice_partitions --help
    usage: simple_seaice_partitions [-h] -m MESHFILENAME -p OUTPUTPREFIX -n
                                    [NPROCSARRAY ...] [-d DATADIR]

    Create sea-ice partitions on LCRC.

    options:
      -h, --help            show this help message and exit
      -m MESHFILENAME, --mesh MESHFILENAME
                            MPAS-Seaice mesh file.
      -p OUTPUTPREFIX, --prefix OUTPUTPREFIX
                            prefix for output partition filenames.
      -n [NPROCSARRAY ...], --nprocs [NPROCSARRAY ...]
                            list of the number of processors to create partition
                            for.
      -d DATADIR, --datadir DATADIR
                            Directory with seaice_QU60km_polar.nc and
                            icePresent_QU60km_polar.nc.

The mesh file is any file that contains the MPAS-Seaice mesh.  Some meshes
are available in `inputdata/share/meshes/mpas/sea-ice` and also each
MPAS-Seaice initial condition in `inputdata/ice/mpas-seaice/<mesh_name>`
contains the MPAS mesh.  Which specific initial condition you choose should
not matter because the mesh should be identical.

The output prefix can be an absolute or relative path prefix for the graph
partition file to be created.  Typically, this will be something like
``partitions/mpas-seaice.graph.info.230313``.  It should end in a date that
matches other existing partition files (i.e. it can't typically be today's
date or E3SM won't find the new partition file) and should not contain the
``.part.<task_count>`` suffix.

You can provide several task counts with ``-n`` for efficiency.  There is a
significant overhead in calling the tool multiple times for different task
counts.

Here is an example:

.. code-block:: bash

    cd /lcrc/group/e3sm/public_html/inputdata/ice/mpas-seaice/WC14to60E2r3
    simple_seaice_partitions -m seaice.WC14to60E2r3.210210.nc -p partitions/mpas-seaice.graph.info.230313 -n 468


Graph partition function
------------------------

A helper function :py:func:`mpas_tools.seaice.partition.gen_seaice_mesh_partition()`
is used within ``create_seaice_partitions``.  It can also be called directly
but must already have the files resulting from ``prepare_seaice_partitions``
available in the output directory.
