
import numpy as np
import os
import xarray
import argparse
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee


def sort_fwd(data, fwd):
    vals = data.values
    vals = vals[fwd - 1]
    return vals


def sort_rev(data, rev):
    vals = data.values
    mask = vals > 0
    vals[mask] = rev[vals[mask] - 1]
    return vals


def sort_mesh(mesh):
    """
    SORT-MESH: sort cells, edges and duals in the mesh
    to improve cache-locality.

    """
    # Authors: Darren Engwirda

    ncells = mesh.dims["nCells"]
    nedges = mesh.dims["nEdges"]
    nduals = mesh.dims["nVertices"]

    cell_fwd = np.arange(0, ncells) + 1
    cell_rev = np.arange(0, ncells) + 1
    edge_fwd = np.arange(0, nedges) + 1
    edge_rev = np.arange(0, nedges) + 1
    dual_fwd = np.arange(0, nduals) + 1
    dual_rev = np.arange(0, nduals) + 1

    # sort cells via RCM ordering of adjacency matrix

    cell_fwd = reverse_cuthill_mckee(cell_del2(mesh)) + 1

    cell_rev = np.zeros(ncells, dtype=np.int32)
    cell_rev[cell_fwd - 1] = np.arange(ncells) + 1

    mesh["cellsOnCell"][:] = \
        sort_rev(mesh["cellsOnCell"], cell_rev)
    mesh["cellsOnEdge"][:] = \
        sort_rev(mesh["cellsOnEdge"], cell_rev)
    mesh["cellsOnVertex"][:] = \
        sort_rev(mesh["cellsOnVertex"], cell_rev)

    for var in mesh.keys():
        dims = mesh.variables[var].dims
        if ("nCells" in dims):
            mesh[var][:] = sort_fwd(mesh[var], cell_fwd)

    mesh["indexToCellID"][:] = np.arange(ncells) + 1

    # sort duals via pseudo-linear cell-wise ordering

    dual_fwd = np.ravel(mesh["verticesOnCell"].values)
    dual_fwd = dual_fwd[dual_fwd > 0]

    __, imap = np.unique(dual_fwd, return_index=True)

    dual_fwd = dual_fwd[np.sort(imap)]

    dual_rev = np.zeros(nduals, dtype=np.int32)
    dual_rev[dual_fwd - 1] = np.arange(nduals) + 1

    mesh["verticesOnCell"][:] = \
        sort_rev(mesh["verticesOnCell"], dual_rev)

    mesh["verticesOnEdge"][:] = \
        sort_rev(mesh["verticesOnEdge"], dual_rev)

    for var in mesh.keys():
        dims = mesh.variables[var].dims
        if ("nVertices" in dims):
            mesh[var][:] = sort_fwd(mesh[var], dual_fwd)

    mesh["indexToVertexID"][:] = np.arange(nduals) + 1

    # sort edges via pseudo-linear cell-wise ordering

    edge_fwd = np.ravel(mesh["edgesOnCell"].values)
    edge_fwd = edge_fwd[edge_fwd > 0]

    __, imap = np.unique(edge_fwd, return_index=True)

    edge_fwd = edge_fwd[np.sort(imap)]

    edge_rev = np.zeros(nedges, dtype=np.int32)
    edge_rev[edge_fwd - 1] = np.arange(nedges) + 1

    mesh["edgesOnCell"][:] = \
        sort_rev(mesh["edgesOnCell"], edge_rev)

    mesh["edgesOnEdge"][:] = \
        sort_rev(mesh["edgesOnEdge"], edge_rev)

    mesh["edgesOnVertex"][:] = \
        sort_rev(mesh["edgesOnVertex"], edge_rev)

    for var in mesh.keys():
        dims = mesh.variables[var].dims
        if ("nEdges" in dims):
            mesh[var][:] = sort_fwd(mesh[var], edge_fwd)

    mesh["indexToEdgeID"][:] = np.arange(nedges) + 1

    return mesh


def cell_del2(mesh):
    """
    CELL-DEL2: form cell-to-cell sparse adjacency graph

    """
    xvec = np.array([], dtype=np.int8)
    ivec = np.array([], dtype=np.int32)
    jvec = np.array([], dtype=np.int32)

    topolOnCell = mesh["nEdgesOnCell"].values
    cellsOnCell = mesh["cellsOnCell"].values

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


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--mesh-file", dest="mesh_file", type=str,
        required=True, help="Path+name to unsorted mesh file.")

    parser.add_argument(
        "--sort-file", dest="sort_file", type=str,
        required=True, help="Path+name to sorted output file.")

    args = parser.parse_args()

    mesh = xarray.open_dataset(args.mesh_file)

    sort_mesh(mesh)

    with open(os.path.join(os.path.dirname(
              args.sort_file), "graph.info"), "w") as fptr:
        cellsOnCell = mesh["cellsOnCell"].values

        ncells = mesh.dims["nCells"]
        nedges = np.count_nonzero(cellsOnCell) // 2

        fptr.write("{} {}\n".format(ncells, nedges))
        for cell in range(ncells):
            data = cellsOnCell[cell, :]
            data = data[data > 0]
            for item in data:
                fptr.write("{} ".format(item))
            fptr.write("\n")

    mesh.to_netcdf(args.sort_file, format="NETCDF4")
