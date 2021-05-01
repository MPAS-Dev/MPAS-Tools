import multiprocessing
import numpy as np
import subprocess
from netCDF4 import Dataset
import argparse
from collections import defaultdict


def create_pool(process_count=None, method='forkserver'):
    """
    Crate a pool for creating masks with Python multiprocessing.  This should
    be called only once at the beginning of the script performing cell culling.
    ``pool.terminate()`` should be called before exiting the script.

    Parameters
    ----------
    process_count : int, optional
        The number of processors or None to use all available processors

    method : {'fork', 'spawn', 'forkserver'}
        The mutiprocessing method

    Returns
    -------
    pool : multiprocessing.Pool
        A pool to use for python-based mask creation.
    """
    pool = None
    multiprocessing.set_start_method(method)
    if process_count is None:
        process_count = multiprocessing.cpu_count()
    else:
        process_count = min(process_count, multiprocessing.cpu_count())

    if process_count > 1:
        pool = multiprocessing.Pool(process_count)

    return pool


def make_partition_files(mesh_filename, num_blocks, num_procs, metis='gpmetis',
                         weight_field=None):
    """
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

    Parameters
    ----------
    mesh_filename : str
        An MPAS mesh file from which the partition should be created

    num_blocks : int
        The number of blocks to decompose

    num_procs : int
        The number of processor to decompose

    metis : str, optional
        The name or path of the Metis executable

    weight_field: str, optional
        The variable in ``mesh_filename`` by which the block decomposition
        should be weighted
    """
    if weight_field is None:
        print("Weight field missing. Defaulting to unweighted graphs.")
        weighted_parts = False
    else:
        weighted_parts = True

    if num_procs > 1:
        proc_decomp = True
    else:
        proc_decomp = False

    with Dataset(mesh_filename, 'r') as grid:
        nCells = len(grid.dimensions['nCells'])

        nEdgesOnCell = grid.variables['nEdgesOnCell'][:]
        cellsOnCell = grid.variables['cellsOnCell'][:] - 1
        if weighted_parts:
            weights = grid.variables[weight_field][:]

    local_num_blocks = 0
    owning_block = [0] * nCells
    block_owner = [0] * int(num_blocks)

    nEdges = 0
    for i in np.arange(0, nCells):
        for j in np.arange(0, nEdgesOnCell[i]):
            if cellsOnCell[i][j] != -1:
                nEdges = nEdges + 1

    nEdges = nEdges // 2

    graph = open('graph.info', 'w+')
    graph.write('%s %s\n' % (nCells, nEdges))
    if weighted_parts:
        wgraph = open('weighted.graph.info', 'w+')
        wgraph.write('%s %s 010\n' % (nCells, nEdges))
    else:
        wgraph = None

    for i in np.arange(0, nCells):
        if weighted_parts:
            wgraph.write('%s ' % int(weights[i]))

        for j in np.arange(0, nEdgesOnCell[i]):
            if weighted_parts:
                if cellsOnCell[i][j] >= 0:
                    wgraph.write('%s ' % (cellsOnCell[i][j] + 1))

            if cellsOnCell[i][j] >= 0:
                graph.write('%s ' % (cellsOnCell[i][j] + 1))
        graph.write('\n')

        if weighted_parts:
            wgraph.write('\n')
    graph.close()

    if weighted_parts:
        wgraph.close()

    command = metis
    if weighted_parts:
        arg1 = "weighted.graph.info"
    else:
        arg1 = "graph.info"
    arg2 = "{}".format(num_blocks)
    subprocess.check_call([command, arg1, arg2], stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL)

    if weighted_parts:
        graph = open('weighted.graph.info.part.{}'.format(num_blocks), 'r')
    else:
        graph = open('graph.info.part.{}'.format(num_blocks), 'r')
    i = -1
    for block in iter(lambda: graph.readline(), ""):
        if i >= 0:
            block_arr = block.split()
            owning_block[i] = int(block_arr[0]) + 1
            local_num_blocks = max(local_num_blocks, owning_block[i])

        i = i + 1
    graph.close()

    if proc_decomp:
        blocksOnBlock = defaultdict(list)
        nEdges = 0

        for i in np.arange(0, nCells):
            for j in np.arange(0, nEdgesOnCell[i]):
                iCell = cellsOnCell[i][j]
                try:
                    can_add = True
                    for block in blocksOnBlock[owning_block[i]]:
                        if block == owning_block[iCell]:
                            can_add = False
                except:
                    can_add = True

                if iCell == -1:
                    can_add = False

                if owning_block[iCell] == owning_block[i]:
                    can_add = False

                if owning_block[iCell] <= 0:
                    can_add = False

                if owning_block[i] <= 0:
                    can_add = False

                if can_add:
                    nEdges = nEdges + 1
                    blocksOnBlock[owning_block[i]].append(owning_block[iCell])

        del blocksOnBlock[0]

        block_graph = open('block.graph.info', 'w+')
        block_graph.write('{} {}\n'.format(local_num_blocks, nEdges // 2))
        for i in np.arange(1, local_num_blocks + 1):
            for block in blocksOnBlock[i]:
                block_graph.write('{} '.format(block))
            block_graph.write('\n')

        block_graph.close()

        command = metis
        arg1 = "block.graph.info"
        arg2 = "{}".format(num_procs)
        subprocess.check_call([command, arg1, arg2], stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)

        block_graph = open('block.graph.info.part.{}'.format(num_procs), 'r')
        iBlock = 0
        for block in iter(lambda: block_graph.readline(), ""):
            block_arr = block.split()
            block_owner[iBlock] = int(block_arr[0])
            iBlock = iBlock + 1

        block_graph.close()

        block_location = open('int_ext_blocks.dat', 'w+')
        interior_blocks = 0
        exterior_blocks = 0
        for i in np.arange(1, local_num_blocks + 1):
            owner = block_owner[i - 1]
            interior = True
            for block in blocksOnBlock[i - 1]:
                if block_owner[block - 1] != owner:
                    interior = False

            if interior:
                block_location.write('1\n')
                interior_blocks = interior_blocks + 1
            else:
                block_location.write('0\n')
                exterior_blocks = exterior_blocks + 1

        block_location.close()

        print('Interior blocks: ', interior_blocks)
        print('Exterior blocks: ', exterior_blocks)


def make_partition_files_entry_point():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--file", dest="mesh_filename", required=True,
                        help="Path to grid file", metavar="FILE")
    parser.add_argument("-m", "--metis", dest="metis", default='gpmetis',
                        help="Path or name of metis executable, default is "
                             "'gpmetis'",
                        metavar="METIS")
    parser.add_argument("-p", "--procs", dest="num_procs", type=int,
                        required=True,
                        help="Number of processors for decomposition",
                        metavar="PROCS")
    parser.add_argument("-b", "--blocks", dest="num_blocks", type=int,
                        required=True,
                        help="Number of blocks for decomposition",
                        metavar="BLOCKS")
    parser.add_argument("-w", "--weights", dest="weight_field",
                        help="Field to weight block partition file on",
                        metavar="VAR")

    args = parser.parse_args()

    make_partition_files(mesh_filename=args.mesh_filename,
                         num_blocks=args.num_blocks, num_procs=args.num_procs,
                         metis=args.metis, weight_field=args.weight_field)
