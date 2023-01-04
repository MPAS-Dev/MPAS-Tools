.. _seaice_partitions:

Graph partitioning
==================

The :py:mod:`mpas_tools.seaice.partition` module contains a function for
creating graph partition files for MPAS-Seaice that are better balanced than
those created from Metis tools directly.

.. _seaice_partitions_:


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

Graph partition function
------------------------

A helper function :py:func:`mpas_tools.seaice.partition.gen_seaice_mesh_partition()`
is used within ``create_seaice_partitions``.  It can also be called directly
but must already have the files resulting from ``prepare_seaice_partitions``
available in the output directory.
