_mesh_vars = [
    'areaCell',
    'cellsOnCell',
    'edgesOnCell',
    'indexToCellID',
    'latCell',
    'lonCell',
    'nEdgesOnCell',
    'verticesOnCell',
    'xCell',
    'yCell',
    'zCell',
    'angleEdge',
    'cellsOnEdge',
    'dcEdge',
    'dvEdge',
    'edgesOnEdge',
    'indexToEdgeID',
    'latEdge',
    'lonEdge',
    'nEdgesOnCell',
    'nEdgesOnEdge',
    'verticesOnEdge',
    'xEdge',
    'yEdge',
    'zEdge',
    'areaTriangle',
    'cellsOnVertex',
    'edgesOnVertex',
    'indexToVertexID',
    'kiteAreasOnVertex',
    'latVertex',
    'lonVertex',
    'xVertex',
    'yVertex',
    'zVertex',
    'weightsOnEdge',
]


def _get_ds_mesh(ds):
    """
    Extract the mesh variables from an xarray Dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        The xarray Dataset containing the mesh data.

    Returns
    -------
    xarray.Dataset
        A new Dataset containing only the specified mesh variables.
    """
    ds_mesh = ds[_mesh_vars]
    return ds_mesh
