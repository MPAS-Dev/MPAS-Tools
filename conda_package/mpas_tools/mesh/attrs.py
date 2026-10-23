"""
CF metadata for MPAS mesh variables and files

``MESH_VAR_ATTRS`` holds the ``long_name``, ``units`` and ``standard_name``
attributes of the MPAS mesh variables.  Units are in the plain form that
udunits parses and CF uses (``m``, ``m2``, ``radians``, ``1`` for
dimensionless).  Index and connectivity variables and masks have no units.

``MpasMeshConverter.x`` and ``MpasCellCuller.x`` write the same attributes,
so the two need to be kept in sync.
"""

CF_VERSION = 'CF-1.8'

MESH_VAR_ATTRS = {
    # global indices
    'indexToCellID': {'long_name': 'global index of each cell'},
    'indexToEdgeID': {'long_name': 'global index of each edge'},
    'indexToVertexID': {'long_name': 'global index of each vertex'},
    # connectivity
    'cellsOnCell': {'long_name': 'cells that neighbor each cell'},
    'nEdgesOnCell': {'long_name': 'number of edges that border each cell'},
    'edgesOnCell': {'long_name': 'edges that border each cell'},
    'verticesOnCell': {'long_name': 'vertices that border each cell'},
    'cellsOnEdge': {'long_name': 'cells that straddle each edge'},
    'edgesOnEdge': {
        'long_name': 'edges that border the cells that straddle each edge'
    },
    'nEdgesOnEdge': {
        'long_name': 'number of edges that border the cells that straddle '
        'each edge'
    },
    'verticesOnEdge': {'long_name': 'vertices that straddle each edge'},
    'cellsOnVertex': {'long_name': 'cells that share each vertex'},
    'edgesOnVertex': {'long_name': 'edges that share each vertex'},
    # coordinates
    'xCell': {'long_name': 'x coordinate of cell centers', 'units': 'm'},
    'yCell': {'long_name': 'y coordinate of cell centers', 'units': 'm'},
    'zCell': {'long_name': 'z coordinate of cell centers', 'units': 'm'},
    'latCell': {
        'long_name': 'latitude of cell centers',
        'units': 'radians',
        'standard_name': 'latitude',
    },
    'lonCell': {
        'long_name': 'longitude of cell centers',
        'units': 'radians',
        'standard_name': 'longitude',
    },
    'xEdge': {'long_name': 'x coordinate of edge midpoints', 'units': 'm'},
    'yEdge': {'long_name': 'y coordinate of edge midpoints', 'units': 'm'},
    'zEdge': {'long_name': 'z coordinate of edge midpoints', 'units': 'm'},
    'latEdge': {
        'long_name': 'latitude of edge midpoints',
        'units': 'radians',
        'standard_name': 'latitude',
    },
    'lonEdge': {
        'long_name': 'longitude of edge midpoints',
        'units': 'radians',
        'standard_name': 'longitude',
    },
    'xVertex': {'long_name': 'x coordinate of vertices', 'units': 'm'},
    'yVertex': {'long_name': 'y coordinate of vertices', 'units': 'm'},
    'zVertex': {'long_name': 'z coordinate of vertices', 'units': 'm'},
    'latVertex': {
        'long_name': 'latitude of vertices',
        'units': 'radians',
        'standard_name': 'latitude',
    },
    'lonVertex': {
        'long_name': 'longitude of vertices',
        'units': 'radians',
        'standard_name': 'longitude',
    },
    # geometry
    'areaCell': {
        'long_name': 'area of each cell in the primal mesh',
        'units': 'm2',
        'standard_name': 'cell_area',
    },
    'areaTriangle': {
        'long_name': 'area of each triangle in the dual mesh',
        'units': 'm2',
    },
    'kiteAreasOnVertex': {
        'long_name': 'area of the part of each dual cell in each cell on the '
        'vertex',
        'units': 'm2',
    },
    'dvEdge': {
        'long_name': 'distance between the vertices at the ends of each edge',
        'units': 'm',
    },
    'dcEdge': {
        'long_name': 'distance between the centers of the cells on each edge',
        'units': 'm',
    },
    'angleEdge': {
        'long_name': 'angle between the normal of each edge and local east',
        'units': 'radians',
    },
    'weightsOnEdge': {
        'long_name': 'weights for reconstructing tangential velocity from '
        'edges on edge',
        'units': '1',
    },
    'meshDensity': {
        'long_name': 'value of the density function used to generate the mesh',
        'units': '1',
    },
    # Coriolis
    'fCell': {
        'long_name': 'Coriolis parameter at cell centers',
        'units': 'radians s-1',
    },
    'fEdge': {
        'long_name': 'Coriolis parameter at edges',
        'units': 'radians s-1',
    },
    'fVertex': {
        'long_name': 'Coriolis parameter at vertices',
        'units': 'radians s-1',
    },
    # masks
    'cullCell': {
        'long_name': 'mask of cells to be removed by the cell culler',
    },
    'boundaryVertex': {
        'long_name': 'mask of vertices with at least one inactive neighboring '
        'cell',
    },
    'boundaryEdge': {
        'long_name': 'mask of edges with only one active neighboring cell',
    },
    'boundaryCell': {
        'long_name': 'mask of cells with at least one inactive neighboring '
        'cell',
    },
    # mesh quality from the mesh converter
    'cellQuality': {
        'long_name': 'ratio of the shortest to the longest edge of each cell',
        'units': '1',
    },
    'triangleQuality': {
        'long_name': 'ratio of the shortest to the longest edge of each dual '
        'triangle',
        'units': '1',
    },
    'triangleAngleQuality': {
        'long_name': 'ratio of the smallest to the largest angle of each dual '
        'triangle',
        'units': '1',
    },
    'obtuseTriangle': {
        'long_name': 'mask of dual triangles with an obtuse angle',
    },
    'gridSpacing': {
        'long_name': 'mean distance from each cell center to its neighbors',
        'units': 'm',
    },
}


def cf_conventions(conventions=None):
    """
    Get a ``Conventions`` attribute for an MPAS mesh file that includes CF

    Entries in ``conventions`` (including any CF version) are kept, and
    ``CF-1.8`` and ``MPAS`` are added if they are missing.

    Parameters
    ----------
    conventions : str, optional
        The existing ``Conventions`` attribute, with entries separated by
        blanks or commas

    Returns
    -------
    conventions : str
        The ``Conventions`` attribute with a CF entry first
    """
    if conventions is None:
        conventions = ''
    entries = conventions.replace(',', ' ').split()
    cf = [entry for entry in entries if entry.startswith('CF-')]
    other = [entry for entry in entries if not entry.startswith('CF-')]
    if 'MPAS' not in other:
        other.append('MPAS')
    cf_entry = cf[0] if cf else CF_VERSION
    return ' '.join([cf_entry] + other)


def add_mesh_attrs(ds):
    """
    Add CF metadata to the mesh variables in a dataset and add CF to its
    ``Conventions``.  Attributes that a variable already has are kept.

    Parameters
    ----------
    ds : xarray.Dataset
        An MPAS mesh dataset, modified in place

    Returns
    -------
    ds : xarray.Dataset
        The same dataset, for convenience
    """
    for name, attrs in MESH_VAR_ATTRS.items():
        if name not in ds:
            continue
        var_attrs = ds[name].attrs
        for key, value in attrs.items():
            var_attrs.setdefault(key, value)
    ds.attrs['Conventions'] = cf_conventions(ds.attrs.get('Conventions'))
    return ds
