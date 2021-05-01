#!/usr/bin/env python

"""
usage: make_partition_files.py [-h] -f FILE [-m METIS] -p PROCS -b BLOCKS
                               [-w VAR]

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Path to grid file
  -m METIS, --metis METIS
                        Path or name of metis executable, default is 'gpmetis'
  -p PROCS, --procs PROCS
                        Number of processors for decomposition
  -b BLOCKS, --blocks BLOCKS
                        Number of blocks for decomposition
  -w VAR, --weights VAR
                        Field to weight block partition file on.

Create weighted and unweighted hierarchical decompositions for use with
multiple blocks per MPI process in any MPAS core.

After running the script, several files are generated which are described
below:

- ``graph.info`` - Single processor graph decomposition file.

- (optional) ``weighted.graph.info`` - ``graph.info`` file used for
  creating weighted partition files.

- ``(weighted.)graph.info.part.BLOCKS`` - partition file that is weighted
  or unweighted and has BLOCKS number of partitions.

- ``block.graph.info`` - graph file of blocks. Equivalent to
  ``blocksOnBlock``

- ``block.graph.info.part.PROCS`` - partition file that has PROCS number of
   partitions.

These can be used in MPAS as the block decomposition file
(``graph.info.part.BLOCKS``) and the proc decomposition file
(``block.graph.info.part.PROCS``) to control the topology of blocks on MPI
tasks.

Additional ``graph.info.part`` files can be created by running the script
multiple times, or by running metis on the graph.info file.

Additional block.graph.info.part files can be created through the same
process, but they can only be used in tandem with the corresponding
``graph.info.part.BLOCKS`` file.
"""

from mpas_tools.parallel import make_partition_files_entry_point


if __name__ == '__main__':
    make_partition_files_entry_point()
