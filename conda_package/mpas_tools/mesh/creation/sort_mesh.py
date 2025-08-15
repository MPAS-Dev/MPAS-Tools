import argparse
import os

import numpy as np
import xarray
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from mpas_tools.io import write_netcdf


def sort_mesh(mesh):
    """
    Sort cells, edges and duals in the mesh
    to improve cache-locality

    Parameters
    ----------
    mesh : xarray.Dataset
        A dataset containing an MPAS mesh to sort

    Returns
    -------
    mesh : xarray.Dataset
        A dataset containing the sorted MPAS mesh
    """
    # Authors: Darren Engwirda

    ncells = mesh.sizes['nCells']
    nedges = mesh.sizes['nEdges']
    nduals = mesh.sizes['nVertices']

    cell_fwd = np.arange(0, ncells) + 1
    cell_rev = np.arange(0, ncells) + 1
    edge_fwd = np.arange(0, nedges) + 1
    edge_rev = np.arange(0, nedges) + 1
    dual_fwd = np.arange(0, nduals) + 1
    dual_rev = np.arange(0, nduals) + 1

    # sort cells via RCM ordering of adjacency matrix

    cell_fwd = reverse_cuthill_mckee(_cell_del2(mesh)) + 1

    cell_rev = np.zeros(ncells, dtype=np.int32)
    cell_rev[cell_fwd - 1] = np.arange(ncells) + 1

    mesh['cellsOnCell'][:] = _sort_rev(mesh['cellsOnCell'], cell_rev)
    mesh['cellsOnEdge'][:] = _sort_rev(mesh['cellsOnEdge'], cell_rev)
    mesh['cellsOnVertex'][:] = _sort_rev(mesh['cellsOnVertex'], cell_rev)

    for var in mesh.keys():
        dims = mesh.variables[var].dims
        if 'nCells' in dims:
            mesh[var][:] = _sort_fwd(mesh[var], cell_fwd)

    mesh['indexToCellID'][:] = np.arange(ncells) + 1

    # sort duals via pseudo-linear cell-wise ordering

    dual_fwd = np.ravel(mesh['verticesOnCell'].values)
    dual_fwd = dual_fwd[dual_fwd > 0]

    __, imap = np.unique(dual_fwd, return_index=True)

    dual_fwd = dual_fwd[np.sort(imap)]

    dual_rev = np.zeros(nduals, dtype=np.int32)
    dual_rev[dual_fwd - 1] = np.arange(nduals) + 1

    mesh['verticesOnCell'][:] = _sort_rev(mesh['verticesOnCell'], dual_rev)

    mesh['verticesOnEdge'][:] = _sort_rev(mesh['verticesOnEdge'], dual_rev)

    for var in mesh.keys():
        dims = mesh.variables[var].dims
        if 'nVertices' in dims:
            mesh[var][:] = _sort_fwd(mesh[var], dual_fwd)

    mesh['indexToVertexID'][:] = np.arange(nduals) + 1

    # sort edges via pseudo-linear cell-wise ordering

    edge_fwd = np.ravel(mesh['edgesOnCell'].values)
    edge_fwd = edge_fwd[edge_fwd > 0]

    __, imap = np.unique(edge_fwd, return_index=True)

    edge_fwd = edge_fwd[np.sort(imap)]

    edge_rev = np.zeros(nedges, dtype=np.int32)
    edge_rev[edge_fwd - 1] = np.arange(nedges) + 1

    mesh['edgesOnCell'][:] = _sort_rev(mesh['edgesOnCell'], edge_rev)

    mesh['edgesOnEdge'][:] = _sort_rev(mesh['edgesOnEdge'], edge_rev)

    mesh['edgesOnVertex'][:] = _sort_rev(mesh['edgesOnVertex'], edge_rev)

    for var in mesh.keys():
        dims = mesh.variables[var].dims
        if 'nEdges' in dims:
            mesh[var][:] = _sort_fwd(mesh[var], edge_fwd)

    mesh['indexToEdgeID'][:] = np.arange(nedges) + 1

    return mesh


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        '--mesh-file',
        dest='mesh_file',
        type=str,
        required=True,
        help='Path+name to unsorted mesh file.',
    )

    parser.add_argument(
        '--sort-file',
        dest='sort_file',
        type=str,
        required=True,
        help='Path+name to sorted output file.',
    )

    parser.add_argument(
        '--format',
        dest='format',
        type=str,
        required=False,
        default='NETCDF4',
        help='Output format for the sorted mesh file.',
    )

    args = parser.parse_args()

    mesh = xarray.open_dataset(args.mesh_file)

    sort_mesh(mesh)

    with open(
        os.path.join(os.path.dirname(args.sort_file), 'graph.info'), 'w'
    ) as fptr:
        cellsOnCell = mesh['cellsOnCell'].values

        ncells = mesh.sizes['nCells']
        nedges = np.count_nonzero(cellsOnCell) // 2

        fptr.write(f'{ncells} {nedges}\n')
        for cell in range(ncells):
            data = cellsOnCell[cell, :]
            data = data[data > 0]
            for item in data:
                fptr.write(f'{item} ')
            fptr.write('\n')

    write_netcdf(mesh, args.sort_file, format=args.format)


def _sort_fwd(data, fwd):
    """
    Apply a forward permutation to a mesh array

    Parameters
    ----------
    data : array-like
        An MPAS mesh array to permute

    fwd : numpy.ndarray
        An array of integers defining the permutation

    Returns
    -------
    data : numpy.ndarray
        The forward permuted MPAS mesh array
    """
    vals = data.values
    vals = vals[fwd - 1]
    return vals


def _sort_rev(data, rev):
    """
    Apply a reverse permutation to a mesh array

    Parameters
    ----------
    data : array-like
        An MPAS mesh array to permute

    rev : numpy.ndarray
        An array of integers defining the permutation

    Returns
    -------
    data : numpy.ndarray
        The reverse permuted MPAS mesh array
    """
    vals = data.values
    mask = vals > 0
    vals[mask] = rev[vals[mask] - 1]
    return vals


def _cell_del2(mesh):
    """
    Form cell-to-cell sparse adjacency graph

    Parameters
    ----------
    mesh : xarray.Dataset
        A dataset containing an MPAS mesh to sort

    Returns
    -------
    del2 : scipy.sparse.csr_matrix
        The cell-to-cell adjacency graph as a sparse matrix
    """
    xvec = np.array([], dtype=np.int8)
    ivec = np.array([], dtype=np.int32)
    jvec = np.array([], dtype=np.int32)

    topolOnCell = mesh['nEdgesOnCell'].values
    cellsOnCell = mesh['cellsOnCell'].values

    for edge in range(np.max(topolOnCell)):
        # cell-to-cell pairs, if edges exist
        mask = topolOnCell > edge
        idx_self = np.argwhere(mask).ravel()
        idx_next = cellsOnCell[mask, edge] - 1

        # cell-to-cell pairs, if cells exist
        mask = idx_next >= 0
        idx_self = idx_self[mask]
        idx_next = idx_next[mask]

        # dummy matrix values, just topol. needed
        val_edge = np.ones(idx_next.size, dtype=np.int8)

        ivec = np.hstack((ivec, idx_self))
        jvec = np.hstack((jvec, idx_next))
        xvec = np.hstack((xvec, val_edge))

    return csr_matrix((xvec, (ivec, jvec)))


if __name__ == '__main__':
    main()
