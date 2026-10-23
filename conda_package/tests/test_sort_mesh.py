import numpy as np
import xarray

from mpas_tools.mesh.conversion import convert, cull
from mpas_tools.mesh.creation.sort_mesh import sort_mesh

from .util import get_test_data_file


def test_sort_mesh_with_isolated_last_cell():
    ds_mesh = xarray.open_dataset(
        get_test_data_file('mesh.QU.1920km.151026.nc')
    )
    ds_mesh = convert(dsIn=ds_mesh)

    # cull every neighbor of the last cell, so that it is still the last cell
    # of the culled mesh but has no neighbors left
    n_cells = ds_mesh.sizes['nCells']
    last = n_cells - 1
    neighbors = ds_mesh.cellsOnCell.values[last]
    neighbors = neighbors[neighbors > 0] - 1
    cull_cell = np.zeros(n_cells, dtype=int)
    cull_cell[neighbors] = 1
    ds_mesh['cullCell'] = ('nCells', cull_cell)
    ds_culled = cull(dsIn=ds_mesh)
    assert (ds_culled.cellsOnCell.values[-1] == 0).all()

    ds_sorted = sort_mesh(ds_culled.copy(deep=True))

    n_culled = ds_culled.sizes['nCells']
    assert ds_sorted.sizes['nCells'] == n_culled
    np.testing.assert_array_equal(
        ds_sorted.indexToCellID.values, np.arange(n_culled) + 1
    )
    # the isolated cell is still somewhere in the mesh, with no neighbors
    isolated = (ds_sorted.cellsOnCell.values == 0).all(axis=1)
    assert isolated.sum() == 1
