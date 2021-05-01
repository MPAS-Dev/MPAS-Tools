import multiprocessing
import argparse
import xarray


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


def make_graph_file(mesh_filename, graph_filename='graph.info',
                    weight_field=None):
    """
    Make a graph file from the MPAS mesh for use in the Metis graph
    partitioning software

    Parameters
    ----------
     mesh_filename : str
        The name of the input MPAS mesh file

    graph_filename : str, optional
        The name of the output graph file

    weight_field : str
        The name of a variable in the MPAS mesh file to use as a field of
        weights
    """

    with xarray.open_dataset(mesh_filename) as ds:

        nCells = ds.sizes['nCells']

        nEdgesOnCell = ds.nEdgesOnCell.values
        cellsOnCell = ds.cellsOnCell.values - 1
        if weight_field is not None:
            if weight_field in ds:
                raise ValueError('weight_field {} not found in {}'.format(
                    weight_field, mesh_filename))
            weights = ds[weight_field].values
        else:
            weights = None

    nEdges = 0
    for i in range(nCells):
        for j in range(nEdgesOnCell[i]):
            if cellsOnCell[i][j] != -1:
                nEdges = nEdges + 1

    nEdges = nEdges/2

    with open(graph_filename, 'w+') as graph:
        if weights is None:
            graph.write('{} {}\n'.format(nCells, nEdges))

            for i in range(nCells):
                for j in range(0, nEdgesOnCell[i]):
                    if cellsOnCell[i][j] >= 0:
                        graph.write('{} '.format(cellsOnCell[i][j]+1))
                graph.write('\n')
        else:
            graph.write('{} {} 010\n'.format(nCells, nEdges))

            for i in range(nCells):
                graph.write('{} '.format(int(weights[i])))
                for j in range(0, nEdgesOnCell[i]):
                    if cellsOnCell[i][j] >= 0:
                        graph.write('{} '.format(cellsOnCell[i][j] + 1))
                graph.write('\n')


def make_graph_file_entry_point():

    parser = argparse.ArgumentParser(
        description='Make a graph file from the MPAS mesh for use in the Metis'
                    ' graph partitioning software',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--mesh_filename", dest="mesh_filename",
                        help="Path to mesh file", metavar="FILE",
                        required=True)
    parser.add_argument("-o", "--graph_filename", dest="graph_filename",
                        help="File name of the graph file (default "
                             "is graph.info)",
                        metavar="FILE", default='graph.info')
    parser.add_argument("-w", "--weights", dest="weight_field",
                        help="Field in the mesh file to use for weights",
                        metavar="VAR")

    args = parser.parse_args()

    make_graph_file(mesh_filename=args.mesh_filename,
                    graph_filename=args.graph_filename,
                    weight_field=args.weight_field)
