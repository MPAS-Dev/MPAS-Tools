#!/usr/bin/env python3
"""
Correct MALI muFriction using a flux-gate-based discharge mismatch correction.

For each MALI flux gate (see mesh_tools_li/make_flux_gate_masks.py), this script
compares the mean modeled discharge (from a MALI fluxGatesOutput file) against the
matched GEUS (Mankoff et al., 2020) observed discharge over the same calendar period,
and derives a multiplicative correction

    muScale_g = (Fm_g / Fo_g) ** exponent

for each gate. The correction is computed independently per gate (not from regional
totals), seeded onto the MALI cells adjacent to that gate's edges, averaged in log
space where gates share a cell, and then creep-filled across the *entire* connected
mesh (including cells outside the current ice mask -- muFriction there is already an
extrapolated/background field) using mesh cell-neighbor connectivity. Cells not
connected to any gate get muScale = 1. The result multiplies the input muFriction:

    muFriction_new = muFriction_old * muScale

This is an empirical, iterative correction, not an exact inversion of the basal
friction law. Gate matching (MALI gate -> Mankoff/GEUS gate) reuses the same
location-based matching as output_processing_li/plot_flux_gates.py (see
mesh_tools_li/flux_gate_utils.py), and creep-fill reuses
mesh_tools_li/extrapolate_variable.py's extrapolate_into_mask() with its 'mean' method.

Example usage:
    ./adjust_mu_friction_from_flux_gates.py \\
        --mali-flux-file fluxGates.nc \\
        --mesh-file mesh.nc \\
        --mu-file optimized_friction.nc \\
        --flux-gate-mask-file fluxGateMasks.nc \\
        --gpkg gates.gpkg \\
        --mankoff-discharge-file gate_D.csv \\
        --mankoff-metadata-file gate_meta.csv \\
        --start-index 20 \\
        --end-index 39 \\
        --output-file corrected_friction.nc \\
        --diagnostics-file flux_gate_mu_correction.csv

Created: 2026-10-01
@author: Trevor Hillebrand
"""

import argparse
import os
import shutil
import sys

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from extrapolate_variable import extrapolate_into_mask
from flux_gate_utils import (
    load_mali_flux,
    load_mask_gate_ids,
    load_mask_gate_locations,
    match_gates_by_location,
)

SEED_FILL_VALUE = -9999.0


def decimal_year_to_timestamp(decimal_year):
    """Inverse of flux_gate_utils.xtime_to_decimal_year's fixed 365-day convention."""
    year = int(np.floor(decimal_year))
    day_of_year = (decimal_year - year) * 365.0
    return pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=day_of_year)


def compute_observed_mean_flux(gate_D, obs_gate_id, start_ts, end_ts):
    """Mean observed discharge for one gate over [start_ts, end_ts], or None if unusable."""
    if obs_gate_id not in gate_D.columns:
        return None
    window = gate_D.loc[start_ts:end_ts, obs_gate_id]
    if window.empty or window.isna().all():
        return None
    return float(window.mean())


def gate_seed_cells(edge_mask_column, cells_on_edge):
    """Unique 0-based cells adjacent to the edges where edge_mask_column == 1."""
    edge_indices = np.flatnonzero(edge_mask_column == 1)
    cells = cells_on_edge[edge_indices, :].ravel() - 1
    cells = np.unique(cells[cells >= 0])
    return edge_indices, cells


def bfs_fill(seed_mask, cell_neighbors, neighbor_valid):
    """
    Multi-source BFS over mesh cell connectivity, starting from seed_mask.

    Returns
    -------
    visited : bool ndarray (nCells,), True for cells reachable from any seed cell
    n_iterations : int, number of BFS rings that filled at least one new cell
    """
    visited = seed_mask.copy()
    frontier = np.flatnonzero(seed_mask)
    n_iterations = 0
    while frontier.size:
        candidates = np.unique(cell_neighbors[frontier][neighbor_valid[frontier]])
        candidates = candidates[~visited[candidates]]
        if candidates.size == 0:
            break
        n_iterations += 1
        visited[candidates] = True
        frontier = candidates
    return visited, n_iterations


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--mali-flux-file', required=True,
                        help='MALI fluxGatesOutput file with iceFluxThroughGates and xtime')
    parser.add_argument('--mesh-file', required=True,
                        help='MALI mesh file with cellsOnCell, nEdgesOnCell, cellsOnEdge, '
                             'xCell, yCell')
    parser.add_argument('--mu-file', required=True,
                        help='File containing muFriction to correct (may be the same as '
                             '--mesh-file, but is not required to be)')
    parser.add_argument('--flux-gate-mask-file', required=True,
                        help='Flux-gate mask file with fluxGateIds and fluxGateEdgeMasks '
                             '(from make_flux_gate_masks.py)')
    parser.add_argument('--gpkg', required=True,
                        help='GeoPackage used to build the mask file, for matching mask '
                             'gates to --mankoff-metadata-file by location')
    parser.add_argument('--match-tolerance', type=float, default=1000.0,
                        help='Max distance (m) between mask and Mankoff gate centroids to '
                             'consider them the same physical gate (default: 1000)')
    parser.add_argument('--mankoff-discharge-file', required=True,
                        help='GEUS gate_D.csv (discharge by gate, Gt/yr)')
    parser.add_argument('--mankoff-metadata-file', required=True,
                        help='GEUS gate_meta.csv (gate locations, regions, names)')
    parser.add_argument('--start-index', type=int, required=True,
                        help='First MALI time index (inclusive) to average modeled flux over')
    parser.add_argument('--end-index', type=int, required=True,
                        help='Last MALI time index (inclusive) to average modeled flux over')
    parser.add_argument('--exponent', type=float, default=1.0 / 3.0,
                        help='Exponent applied to Fm/Fo to get muScale (default: 1/3)')
    parser.add_argument('--min-scale', type=float, default=None,
                        help='Lower bound applied to the final muScale field (default: none)')
    parser.add_argument('--max-scale', type=float, default=None,
                        help='Upper bound applied to the final muScale field (default: none)')
    parser.add_argument('--output-file', required=True,
                        help='Output NetCDF file with corrected muFriction (a copy of '
                             '--mu-file with muFriction and diagnostic fields updated)')
    parser.add_argument('--diagnostics-file', default=None,
                        help='Optional output CSV with a gate-by-gate diagnostic table')
    args = parser.parse_args(argv)

    for path in [args.mali_flux_file, args.mesh_file, args.mu_file,
                args.flux_gate_mask_file, args.gpkg, args.mankoff_discharge_file,
                args.mankoff_metadata_file]:
        if not os.path.exists(path):
            parser.error(f"File not found: {path}")
    if os.path.exists(args.output_file):
        parser.error(f"Output file exists: {args.output_file}")

    print(f"Reading gate IDs from mask file: {args.flux_gate_mask_file}")
    gate_ids = load_mask_gate_ids(args.flux_gate_mask_file)

    print(f"Reading MALI flux output: {args.mali_flux_file}")
    years, flux = load_mali_flux([args.mali_flux_file])
    if flux.shape[1] != len(gate_ids):
        sys.exit(f"ERROR: nFluxGates mismatch between flux file ({flux.shape[1]}) "
                 f"and mask file ({len(gate_ids)})")
    n_time = flux.shape[0]
    if not (0 <= args.start_index <= args.end_index < n_time):
        sys.exit(f"ERROR: require 0 <= --start-index <= --end-index < nTime ({n_time}); "
                 f"got start-index={args.start_index}, end-index={args.end_index}")

    Fm = flux[args.start_index:args.end_index + 1, :].mean(axis=0)
    start_ts = decimal_year_to_timestamp(years[args.start_index])
    end_ts = decimal_year_to_timestamp(years[args.end_index])
    print(f"Averaging modeled flux over indices [{args.start_index}, {args.end_index}] "
         f"({start_ts.date()} to {end_ts.date()})")

    print(f"Reading Mankoff metadata/discharge: {args.mankoff_metadata_file}, "
         f"{args.mankoff_discharge_file}")
    meta = pd.read_csv(args.mankoff_metadata_file, index_col='gate')
    gate_D = pd.read_csv(args.mankoff_discharge_file, index_col='Date', parse_dates=True)
    gate_D.columns = gate_D.columns.astype(int)
    if gate_D.loc[start_ts:end_ts].empty:
        sys.exit(f"ERROR: requested period {start_ts.date()} to {end_ts.date()} does not "
                 f"overlap Mankoff discharge data ({gate_D.index.min().date()} to "
                 f"{gate_D.index.max().date()})")

    mask_locations = load_mask_gate_locations(args.gpkg, gate_ids)
    match = match_gates_by_location(mask_locations, meta, tolerance=args.match_tolerance)

    print(f"Reading mesh connectivity: {args.mesh_file}")
    with Dataset(args.mesh_file, 'r') as f:
        n_cells = len(f.dimensions['nCells'])
        cells_on_cell = f.variables['cellsOnCell'][:]
        cells_on_edge = f.variables['cellsOnEdge'][:]
        x_cell = f.variables['xCell'][:]
        y_cell = f.variables['yCell'][:]
    cell_neighbors = cells_on_cell - 1
    neighbor_valid = cell_neighbors >= 0

    with Dataset(args.flux_gate_mask_file, 'r') as f:
        edge_masks = f.variables['fluxGateEdgeMasks'][:, :]
        if edge_masks.shape[0] != cells_on_edge.shape[0]:
            sys.exit(f"ERROR: nEdges mismatch between mask file ({edge_masks.shape[0]}) "
                     f"and mesh file ({cells_on_edge.shape[0]})")

    # Per-gate diagnostics and multiplicative-correction seeding.
    diagnostics = []
    log_mu_scale_sum = np.zeros(n_cells)
    log_mu_scale_count = np.zeros(n_cells, dtype=int)
    for i, mali_gate_id in enumerate(gate_ids):
        edge_indices, seed_cells = gate_seed_cells(edge_masks[:, i], cells_on_edge)
        obs_gate_id = match.get(mali_gate_id)
        region = None
        mankoff_name = None
        Fo = None
        status = None

        if obs_gate_id is None:
            status = 'unmatched (no gate_meta match)'
        else:
            region = meta.loc[obs_gate_id, 'region']
            mankoff_name = meta.loc[obs_gate_id, 'Mouginot_2019']
            Fo = compute_observed_mean_flux(gate_D, obs_gate_id, start_ts, end_ts)
            if Fo is None:
                status = 'no observation in window'

        Fm_g = Fm[i]
        ratio = np.nan
        mu_scale = np.nan
        log_mu_scale = np.nan
        if status is None:
            if not np.isfinite(Fm_g):
                status = 'Fm non-finite'
            elif Fm_g <= 0:
                status = 'Fm<=0'
            elif not np.isfinite(Fo):
                status = 'Fo non-finite'
            elif Fo <= 0:
                status = 'Fo<=0'
            else:
                ratio = Fm_g / Fo
                log_mu_scale = args.exponent * np.log(ratio)
                mu_scale = np.exp(log_mu_scale)

        if status is None:
            if seed_cells.size == 0:
                status = 'no adjacent cells (gate fully ambiguous)'
            else:
                status = 'used'
                log_mu_scale_sum[seed_cells] += log_mu_scale
                log_mu_scale_count[seed_cells] += 1

        diagnostics.append({
            'mali_gate_index': i,
            'mali_gate_id': mali_gate_id,
            'mankoff_gate_id': obs_gate_id,
            'mankoff_name': mankoff_name,
            'region': region,
            'n_edges': len(edge_indices),
            'n_seed_cells': len(seed_cells),
            'Fm_Gt_per_yr': Fm_g,
            'Fo_Gt_per_yr': Fo,
            'Fm_over_Fo': ratio,
            'muScale': mu_scale,
            'logMuScale': log_mu_scale,
            'status': status,
        })

    diagnostics_df = pd.DataFrame(diagnostics)
    n_matched = diagnostics_df['mankoff_gate_id'].notna().sum()
    n_used = (diagnostics_df['status'] == 'used').sum()

    print(f"\nGate matching/correction summary:")
    print(f"  {len(gate_ids)} MALI gates, {n_matched} matched to Mankoff gates, "
         f"{n_used} used for the correction, {len(gate_ids) - n_used} skipped")
    for status, count in diagnostics_df['status'].value_counts().items():
        print(f"    {count:4d}  {status}")
    used = diagnostics_df[diagnostics_df['status'] == 'used']
    if not used.empty:
        print(f"  Fm/Fo over used gates: min={used['Fm_over_Fo'].min():.3f}, "
             f"mean={used['Fm_over_Fo'].mean():.3f}, median={used['Fm_over_Fo'].median():.3f}, "
             f"max={used['Fm_over_Fo'].max():.3f}")
        print(f"  gate muScale over used gates: min={used['muScale'].min():.3f}, "
             f"mean={used['muScale'].mean():.3f}, median={used['muScale'].median():.3f}, "
             f"max={used['muScale'].max():.3f}")

    # Seed cells: geometric mean (mean in log space) across gates sharing a cell.
    seed_mask = log_mu_scale_count > 0
    log_mu_scale_seed = np.where(seed_mask, np.divide(
        log_mu_scale_sum, np.where(log_mu_scale_count > 0, log_mu_scale_count, 1)), np.nan)
    print(f"  {int(seed_mask.sum())} seed cells from {n_used} used gate(s)")

    # Creep-fill the correction across the full connected mesh, in log space.
    var_value = np.where(seed_mask, log_mu_scale_seed, 0.0)
    log_mu_scale_filled = extrapolate_into_mask(
        var_value, seed_mask, cell_neighbors, neighbor_valid, x_cell, y_cell, 'mean')

    visited, n_iterations = bfs_fill(seed_mask, cell_neighbors, neighbor_valid)
    disconnected = ~visited
    n_disconnected = int(disconnected.sum())
    if n_disconnected > 0:
        print(f"WARNING: {n_disconnected} cell(s) not connected to any flux-gate seed; "
             f"muScale set to 1 there")
    log_mu_scale_filled[disconnected] = 0.0
    print(f"  Creep fill: {n_iterations} iteration(s), "
         f"{int(visited.sum())} cell(s) filled, {n_disconnected} disconnected")

    mu_scale_final = np.exp(log_mu_scale_filled)
    if args.min_scale is not None:
        mu_scale_final = np.maximum(mu_scale_final, args.min_scale)
    if args.max_scale is not None:
        mu_scale_final = np.minimum(mu_scale_final, args.max_scale)
    log_mu_scale_final = np.log(mu_scale_final)
    print(f"  Final muFrictionScale over all nCells: min={mu_scale_final.min():.3f}, "
         f"mean={mu_scale_final.mean():.3f}, median={np.median(mu_scale_final):.3f}, "
         f"max={mu_scale_final.max():.3f}")

    print(f"Writing output: {args.output_file}")
    shutil.copyfile(args.mu_file, args.output_file)
    with Dataset(args.output_file, 'r+') as f:
        if 'muFriction' not in f.variables:
            sys.exit(f"ERROR: {args.mu_file} has no muFriction variable")
        if len(f.dimensions['nCells']) != n_cells:
            sys.exit(f"ERROR: nCells mismatch between --mu-file ({len(f.dimensions['nCells'])}) "
                     f"and --mesh-file ({n_cells})")
        mu_friction_old = f.variables['muFriction'][0, :]
        f.variables['muFriction'][0, :] = mu_friction_old * mu_scale_final

        scale_var = f.createVariable('muFrictionScale', 'f8', ('nCells',))
        scale_var[:] = mu_scale_final
        scale_var.long_name = 'Multiplicative correction applied to muFriction, from flux-gate discharge mismatch'
        scale_var.units = '1'

        log_scale_var = f.createVariable('logMuFrictionScale', 'f8', ('nCells',))
        log_scale_var[:] = log_mu_scale_final
        log_scale_var.long_name = 'log(muFrictionScale)'
        log_scale_var.units = '1'

        seed_var = f.createVariable('fluxGateLogMuFrictionScaleSeed', 'f8', ('nCells',),
                                    fill_value=SEED_FILL_VALUE)
        seed_var[:] = np.where(seed_mask, log_mu_scale_seed, SEED_FILL_VALUE)
        seed_var.long_name = ('log(muScale) seeded directly from flux-gate edges, before '
                              'creep-fill; fill value elsewhere')

    if args.diagnostics_file:
        diagnostics_df.to_csv(args.diagnostics_file, index=False)
        print(f"Wrote {args.diagnostics_file}")

    print(f"Wrote {args.output_file}")


if __name__ == '__main__':
    main()
