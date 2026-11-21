import numpy as np
import pytest
import xarray as xr

from mpas_tools.mesh.conversion import cull
from mpas_tools.ocean.viz.transect import (
    compute_transect,
    interp_mpas_to_transect_nodes,
)
from mpas_tools.ocean.viz.transect.vert import _pool_adjacent_violators
from mpas_tools.planar_hex import make_planar_hex_mesh
from mpas_tools.translate import center

N_VERT_LEVELS = 10


@pytest.fixture(scope='module')
def ds_mesh():
    ds_mesh = make_planar_hex_mesh(
        nx=20, ny=20, dc=10e3, nonperiodic_x=True, nonperiodic_y=True
    )
    ds_mesh = cull(ds_mesh)
    center(ds_mesh)
    return ds_mesh


def test_compute_transect_matches_prism_model(ds_mesh):
    """
    Interfaces at transect nodes are those of triangular prisms whose vertex
    interfaces are averages over the adjacent cells in which they are valid,
    even where the same interface lies at very different depths in
    neighboring cells
    """
    n_cells = ds_mesh.sizes['nCells']
    rng = np.random.default_rng(seed=0)
    min_level = np.zeros(n_cells, dtype=int)
    max_level = rng.integers(2, N_VERT_LEVELS, n_cells)
    bottom_depth = rng.uniform(100.0, 500.0, n_cells)
    ssh = rng.uniform(-50.0, 0.0, n_cells)
    thickness = rng.uniform(0.2, 1.0, (n_cells, N_VERT_LEVELS))
    k = np.arange(N_VERT_LEVELS)
    valid = k <= max_level[:, np.newaxis]
    thickness = np.where(valid, thickness, 0.0)
    thickness *= ((ssh + bottom_depth) / thickness.sum(axis=1))[:, None]

    ds_transect = _compute(
        ds_mesh, thickness, bottom_depth, min_level, max_level
    )

    z_interface = ds_transect.zTransectNode.values[:, 0::2]
    assert np.all(np.diff(z_interface, axis=1) <= 0.0)

    z_ref = _prism_model_interfaces(
        ds_transect, thickness, bottom_depth, min_level, max_level
    )
    relevant = _interfaces_in_adjacent_segments(ds_transect)
    assert np.allclose(z_interface[relevant], z_ref[relevant], atol=1e-8)

    # interfaces are a weighted average of interfaces in the cells used for
    # interpolation, so the seafloor can't be below all their beds
    cells = ds_transect.interpCellIndices.values
    weights = ds_transect.interpHorizCellWeights.values
    deepest_bed = np.where(weights > 0.0, -bottom_depth[cells], np.inf).min(
        axis=1
    )
    z_seafloor = ds_transect.zSeafloor.values
    has_seafloor = np.isfinite(z_seafloor)
    assert np.all(z_seafloor[has_seafloor] >= deepest_bed[has_seafloor] - 1e-8)


def test_compute_transect_flat_ssh_and_seafloor(ds_mesh):
    n_cells = ds_mesh.sizes['nCells']
    thickness, bottom_depth = _flat_thickness(n_cells)
    min_level = np.zeros(n_cells, dtype=int)
    max_level = (N_VERT_LEVELS - 1) * np.ones(n_cells, dtype=int)
    ds_transect = _compute(
        ds_mesh, thickness, bottom_depth, min_level, max_level
    )
    assert np.allclose(ds_transect.ssh, 0.0)
    assert np.allclose(ds_transect.zSeafloor, -bottom_depth[0])
    z_ref = np.zeros(N_VERT_LEVELS + 1)
    z_ref[1:] = -np.cumsum(thickness[0, :])
    assert np.allclose(ds_transect.zTransectNode.values[:, 0::2], z_ref)


def test_interp_to_nodes_is_linear_in_z(ds_mesh):
    """
    Interpolation to nodes is linear in z between the middles of layers, so
    zMid interpolates to the elevation of the nodes, even where layer
    thicknesses vary with depth
    """
    n_cells = ds_mesh.sizes['nCells']
    thickness, bottom_depth = _flat_thickness(n_cells)
    min_level = np.zeros(n_cells, dtype=int)
    max_level = (N_VERT_LEVELS - 1) * np.ones(n_cells, dtype=int)
    ds_transect = _compute(
        ds_mesh, thickness, bottom_depth, min_level, max_level
    )

    z_top = np.zeros((n_cells, N_VERT_LEVELS + 1))
    z_top[:, 1:] = -np.cumsum(thickness, axis=1)
    z_mid = 0.5 * (z_top[:, 0:-1] + z_top[:, 1:])
    da_z_mid = xr.DataArray(z_mid, dims=('nCells', 'nVertLevels'))
    da_z_interface = xr.DataArray(z_top, dims=('nCells', 'nVertLevelsP1'))

    z_node = ds_transect.zTransectNode.values
    # the top and bottom nodes take the values in the top and bottom layers
    z_mid_nodes = interp_mpas_to_transect_nodes(ds_transect, da_z_mid).values
    assert np.allclose(z_mid_nodes[:, 1:-1], z_node[:, 1:-1])

    z_interface_nodes = interp_mpas_to_transect_nodes(
        ds_transect, da_z_interface
    ).values
    assert np.allclose(z_interface_nodes, z_node)


def test_compute_transect_land_in_middle(ds_mesh):
    """
    Land cells in the middle of a transect make invalid segments rather
    than being removed
    """
    n_cells = ds_mesh.sizes['nCells']
    thickness, bottom_depth = _flat_thickness(n_cells)
    min_level = np.zeros(n_cells, dtype=int)
    max_level = (N_VERT_LEVELS - 1) * np.ones(n_cells, dtype=int)
    land = np.abs(ds_mesh.xCell.values) < 15e3
    max_level[land] = -1

    ds_transect = _compute(
        ds_mesh, thickness, bottom_depth, min_level, max_level
    )
    assert ds_transect.sizes['nHorizNodes'] == (
        ds_transect.sizes['nSegments'] + 1
    )
    cell_indices = ds_transect.horizCellIndices.values
    invalid = cell_indices < 0
    assert np.any(invalid)
    assert not invalid[0] and not invalid[-1]
    assert not np.any(land[cell_indices[~invalid]])
    assert not np.any(ds_transect.validCells.values[invalid, :])


def test_pool_adjacent_violators():
    values = np.array([0.0, -200.0, -400.0, -300.0, -400.0])
    expected = np.array([0.0, -200.0, -350.0, -350.0, -400.0])
    assert np.allclose(_pool_adjacent_violators(values), expected)


def _flat_thickness(n_cells):
    # layers thicken with depth
    thickness_1d = 5.0 + 3.0 * np.arange(N_VERT_LEVELS)
    thickness = np.tile(thickness_1d, (n_cells, 1))
    bottom_depth = thickness.sum(axis=1)
    return thickness, bottom_depth


def _compute(ds_mesh, thickness, bottom_depth, min_level, max_level):
    # avoid lining up with triangle sides
    x = xr.DataArray([-80e3, 5e3, 80e3], dims=('nPoints',))
    y = xr.DataArray([-30e3, 2e3, 35e3], dims=('nPoints',))
    return compute_transect(
        x=x,
        y=y,
        ds_horiz_mesh=ds_mesh,
        layer_thickness=xr.DataArray(
            thickness, dims=('nCells', 'nVertLevels')
        ),
        bottom_depth=xr.DataArray(bottom_depth, dims=('nCells',)),
        min_level_cell=xr.DataArray(min_level, dims=('nCells',)),
        max_level_cell=xr.DataArray(max_level, dims=('nCells',)),
        spherical=False,
    )


def _prism_model_interfaces(
    ds_transect, thickness, bottom_depth, min_level, max_level
):
    """
    A straightforward implementation of interfaces at transect nodes for
    comparison
    """
    n_cells = len(bottom_depth)
    z = np.zeros((n_cells, N_VERT_LEVELS + 1))
    z[:, -1] = -bottom_depth
    for k in range(N_VERT_LEVELS - 1, -1, -1):
        valid = np.logical_and(k >= min_level, k <= max_level)
        z[:, k] = z[:, k + 1] + np.where(valid, thickness[:, k], 0.0)
    kp1 = np.arange(N_VERT_LEVELS + 1)
    valid_interface = np.logical_and(
        kp1 >= min_level[:, np.newaxis], kp1 <= max_level[:, np.newaxis] + 1
    )

    cells = ds_transect.interpCellIndices.values
    weights = ds_transect.interpHorizCellWeights.values
    n_nodes = cells.shape[0]
    z_nodes = np.full((n_nodes, N_VERT_LEVELS + 1), np.nan)
    for node in range(n_nodes):
        total = np.zeros(N_VERT_LEVELS + 1)
        for end in range(2):
            end_cells = cells[node, 3 * end : 3 * end + 3]
            end_weights = weights[node, 3 * end : 3 * end + 3]
            side_weight = end_weights.sum()
            if side_weight == 0.0:
                continue
            values = np.full(N_VERT_LEVELS + 1, np.nan)
            for k in range(N_VERT_LEVELS + 1):
                numer = 0.0
                denom = 0.0
                for cell, weight in zip(end_cells, end_weights, strict=True):
                    if cell >= 0 and valid_interface[cell, k]:
                        numer += weight * z[cell, k]
                        denom += weight
                if denom > 0.0:
                    values[k] = numer / denom
            finite = np.isfinite(values)
            values[finite] = _isotonic_non_increasing(values[finite])
            total += side_weight * values
        z_nodes[node, :] = total
    return z_nodes


def _isotonic_non_increasing(values):
    """
    The closest non-increasing sequence in the least-squares sense, from the
    min-max formula
    """
    x = -np.asarray(values)
    n = len(x)
    cumsum = np.concatenate([[0.0], np.cumsum(x)])
    result = np.zeros(n)
    for i in range(n):
        result[i] = max(
            min((cumsum[k + 1] - cumsum[j]) / (k + 1 - j) for k in range(i, n))
            for j in range(i + 1)
        )
    return -result


def _interfaces_in_adjacent_segments(ds_transect):
    valid = ds_transect.validCells.values[:, 0::2]
    n_segments, n_levels = valid.shape
    relevant = np.zeros((n_segments + 1, n_levels + 1), dtype=bool)
    for dseg in range(2):
        for dk in range(2):
            relevant[dseg : dseg + n_segments, dk : dk + n_levels] |= valid
    return relevant
