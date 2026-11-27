#!/usr/bin/env python
"""
Calibrate the damage threshold for MALI's crevasse-depth calving law.

For each Antarctic ice-shelf region (the ``ice_shelf_extended_*`` features in
geometric_features, tag ``ais_ice_shelf_calving_calibration``) choose the
damage threshold that best preserves the present calving front, given the
``damage`` field along that front in a MALI output file.

Calving-front cells are dynamic floating cells that either neighbor ice-free
ocean, or neighbor a non-dynamic floating cell that itself neighbors ice-free
ocean.  A cell calves when its damage >= threshold, so a threshold T is judged
by simulating the retreat it would cause starting from the observed front:
front cells with damage >= T are removed, which exposes their dynamic floating
neighbors, which are removed in turn if damage >= T, and so on.  (This captures
thin low-damage margins that protect a high-damage interior.)  The calibrated
threshold is the smallest T for which the floating-area loss stays within
--max-area-loss (default 5%), plus an optional --growth-margin (default 0).

By default (--criterion percentile) the threshold is instead set to a
percentile (--percentile, default 90) of the damage on the calving-front cells
(--percentile-basis front), on the floating cells behind the front, i.e. not on
the front but within --near-dist-km (default 50 km) of it (basis behind), or the
larger of those two percentiles (basis max, the default, which is the safer
choice).  Floating cells with damage below --min-damage (default 0.05) are left
out of these statistics and of the histograms.  The floating-area loss that
results is still computed and reported.  --criterion area-loss selects the
area-loss method described above.

Output: a netCDF file with ``damageCalvingThreshold(Time, nCells)``, constant
within each region.  Where regions overlap the largest value is used (least
retreat).  Cells outside all regions get the median of the calibrated values.

Requires numpy, scipy, netCDF4, matplotlib, shapely (>=2), pyproj and
geometric_features.  Example:
  calibrate_damage_threshold_for_cd_calving.py -i output.nc -m mesh.nc \\
      -o damageCalvingThreshold.nc
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import shapely

from calving_calibration_common import (
    add_region_arguments, combine_region_values, draw_base_map, estimate_cell_area,
    finish_map, identify_front_cells, load_regions, neighbors_from_vertices,
    read_bed_topography, region_mask, write_calibration_file)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True,
                   help='MALI output file with cellMask, damage, xCell, yCell, '
                        'dcEdge and cellsOnVertex')
    p.add_argument('-m', '--mesh', default=None,
                   help='MALI mesh file providing bedTopography (same cells as input); '
                        'not needed if the input file has bedTopography')
    p.add_argument('-o', '--output', default='damageCalvingThreshold.nc')
    add_region_arguments(p)
    p.add_argument('--max-area-loss', type=float, default=0.05,
                   help='floating-area fraction that must stay strictly below this '
                        '(default 0.05)')
    p.add_argument('--criterion', choices=['area-loss', 'percentile'],
                   default='percentile',
                   help='how to set each region\'s threshold (default percentile)')
    p.add_argument('--percentile', type=float, default=90.0,
                   help='percentile of damage used by --criterion percentile '
                        '(default 90)')
    p.add_argument('--percentile-basis', choices=['front', 'behind', 'max'],
                   default='max',
                   help='cells the percentile is taken over: calving-front cells, '
                        'floating cells behind the front, or the larger of the two '
                        'percentiles (default max)')
    p.add_argument('--min-damage', type=float, default=0.05,
                   help='floating cells with damage below this are excluded from '
                        'the percentile statistics and histograms (default 0.05)')
    p.add_argument('--near-dist-km', type=float, default=50.0,
                   help='how far behind the calving front (km) the cells '
                        'counted as "behind the front" extend (default 50)')
    p.add_argument('--growth-margin', type=float, default=0.0,
                   help='damage added to the minimal threshold to favor minor '
                        'growth (default 0)')
    p.add_argument('--time-index', type=int, default=0)
    p.add_argument('--fig-dir', default='figures',
                   help='directory for per-region figures (default figures)')
    return p.parse_args()


def region_groups(in_reg, front, floating, damage, x, y, near_dist, min_damage):
    """Indices of the region's calving-front cells and of its other floating
    cells within near_dist (m) of the calving front ("behind the front").
    Cells with damage < min_damage are excluded from both groups."""
    from scipy.spatial import cKDTree
    keep = in_reg & floating & (damage >= min_damage)
    fr_all = np.where(in_reg & floating & front)[0]  # distances use the whole front
    fr = np.where(keep & front)[0]
    cand = np.where(keep & ~front)[0]
    if len(fr_all) == 0 or len(cand) == 0:
        return fr, cand[:0]
    dist, _ = cKDTree(np.column_stack([x[fr_all], y[fr_all]])).query(
        np.column_stack([x[cand], y[cand]]))
    return fr, cand[dist <= near_dist]


def plot_histograms(ax, thr, fr, near, damage, near_dist, label, min_damage):
    """Overlaid histograms of damage on the region's calving-front cells and on
    its other floating-ice cells within near_dist (m) of the calving front."""
    bins = np.linspace(0.0, 1.0, 41)
    ax.hist(damage[fr], bins=bins, density=True, alpha=0.55, color='tab:red',
            label=f'calving front (n={len(fr)})')
    ax.hist(damage[near], bins=bins, density=True, alpha=0.55, color='tab:blue',
            label=f'not on front, within {near_dist / 1e3:.0f} km (n={len(near)})')
    ax.axvline(thr, color='k', ls='--', label=f'chosen threshold = {thr:.3f} ({label})')
    ax.set_xlabel('damage')
    ax.set_ylabel('probability density')
    ax.set_title(f'Damage (floating ice, damage >= {min_damage:g}): '
                 'calving front vs. nearby')
    ax.legend(fontsize=8)


def plot_map(fig, ax, geom, thr, loss, removed, in_reg, front, ice, floating, damage, x, y,
             dc):
    """Map of damage on ice cells near a region, with contours at the chosen
    threshold and at 0.01, 0.05 and 0.1 below it.  The color range starts at the
    lowest damage on the region's floating ice."""
    tri, idx, extent = draw_base_map(ax, geom, in_reg, front, ice, floating, x, y, dc)
    if tri is not None:
        d = damage[idx]
        dmin = damage[in_reg & floating].min()
        pc = ax.tripcolor(tri, d, vmin=dmin, vmax=1, cmap='turbo', shading='gouraud')
        fig.colorbar(pc, ax=ax, label='damage', shrink=0.8)
        for off, color, ls, lw in ((0.1, 'white', ':', 1.0), (0.05, 'white', '--', 1.0),
                                   (0.01, 'white', '-', 1.0), (0.0, 'red', '-', 1.5)):
            lev = thr - off
            if d[floating[idx]].min() < lev < d[floating[idx]].max():
                ax.tricontour(tri, d, levels=[lev], colors=color, linestyles=ls,
                              linewidths=lw)
    fi = np.where(front & in_reg)[0]
    ri = np.where(removed)[0]
    ax.scatter(x[ri], y[ri], s=6, marker='s', c='red', alpha=0.6, linewidths=0)
    ax.scatter(x[fi], y[fi], s=10, c='none', edgecolors='white', linewidths=0.6)
    finish_map(ax, geom, extent)
    ax.scatter([], [], s=30, marker='s', c='0.6', label='grounded ice')
    ax.plot([], [], 'r-', lw=1.5, label=f'damage = {thr:.3f} (chosen threshold)')
    for off, ls in ((0.01, '-'), (0.05, '--'), (0.1, ':')):
        ax.plot([], [], color='white', ls=ls, lw=1.0,
                label=f'damage = {thr - off:.3f} (threshold - {off:g})')
    ax.scatter([], [], s=10, c='none', edgecolors='k', label='calving-front cells (outlined)')
    ax.scatter([], [], s=12, marker='s', c='red', alpha=0.6,
               label=f'cells removed (connected to terminus): {100 * loss:.2f}% of floating area')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fontsize=8, ncol=2,
              facecolor='0.35', labelcolor='white')
    ax.set_title(f'Damage map (area removed {100 * loss:.2f}%)')


def plot_region(name, geom, thr, loss, removed, in_reg, front, ice, floating, damage, x, y,
                dc, fig_dir, fr, near, near_dist, label, min_damage):
    """Two-panel figure: damage map with threshold contours (left) and overlaid
    histograms of front vs. nearby damage (right)."""
    fig, (ax_map, ax_hist) = plt.subplots(1, 2, figsize=(17, 8),
                                          gridspec_kw={'width_ratios': [1.25, 1]})
    plot_map(fig, ax_map, geom, thr, loss, removed, in_reg, front, ice, floating, damage,
             x, y, dc)
    plot_histograms(ax_hist, thr, fr, near, damage, near_dist, label, min_damage)
    fig.suptitle(f'{name}: chosen damageCalvingThreshold = {thr:.3f}', fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=130, bbox_inches='tight')
    plt.close(fig)


def simulate_retreat(thr, dyn_float, damage, area, adj, seeds_exposed, region):
    """Remove region dynamic-floating cells with damage >= thr, iteratively
    starting from cells exposed to ocean.  Returns the boolean removed mask."""
    removable = region & dyn_float & (damage >= thr)
    removed = np.zeros_like(region)
    exposed = seeds_exposed.copy()
    frontier = removable & exposed
    while frontier.any():
        removed |= frontier
        # neighbors of newly removed cells become exposed
        newly = np.asarray(adj[frontier].sum(axis=0)).ravel() > 0
        frontier = removable & newly & ~removed
    return removed


def main():
    args = parse_args()
    ds = netCDF4.Dataset(args.input)
    t = args.time_index
    ncells = len(ds.dimensions['nCells'])
    mask = ds['cellMask'][t, :].astype(int)
    damage = np.asarray(ds['damage'][t, :], dtype=float)
    bed = read_bed_topography(ds, args.mesh, t)
    x = np.asarray(ds['xCell'][:])
    y = np.asarray(ds['yCell'][:])
    dc = float(np.median(ds['dcEdge'][:]))
    # hexagonal Voronoi cell area from the median cell spacing (no areaCell in file)
    area = estimate_cell_area(dc, ncells)
    adj = neighbors_from_vertices(np.asarray(ds['cellsOnVertex'][:]), ncells)

    cls = identify_front_cells(mask, bed, adj)
    ice, floating, dyn_float = cls['ice'], cls['floating'], cls['dyn_float']
    front, exposed0 = cls['front'], cls['exposed0']
    print(f'{front.sum()} calving-front cells, {dyn_float.sum()} dynamic floating cells, '
          f'cell area {area[0] / 1e6:.2f} km^2')

    regions = load_regions(args.geometric_features, args.region_glob, args.region_tag,
                           args.crs)
    os.makedirs(args.fig_dir, exist_ok=True)
    result = {}
    hdr = f"{'region':45s} {'front':>6s} {'d_p10':>6s} {'d_p50':>6s} {'d_p90':>6s} " \
          f"{'thr':>6s} {'loss%':>6s}"
    print(hdr)
    for name, geom in regions.items():
        in_reg = region_mask(geom, x, y)
        fl_area = area[in_reg & floating].sum()
        fr = in_reg & front
        fr_idx, near_idx = region_groups(in_reg, front, floating, damage, x, y,
                                         1.e3 * args.near_dist_km, args.min_damage)
        if fl_area == 0 or fr.sum() == 0 or (args.criterion == 'percentile'
                                             and len(fr_idx) == 0):
            print(f'{name:45s} no floating/front cells -> default later')
            result[name] = (in_reg, np.nan)
            continue
        if args.criterion == 'percentile':
            pcts = {'front': np.percentile(damage[fr_idx], args.percentile)}
            if len(near_idx) > 0:
                pcts['behind'] = np.percentile(damage[near_idx], args.percentile)
            if args.percentile_basis == 'max':
                thr_min = max(pcts.values())
            elif args.percentile_basis in pcts:
                thr_min = pcts[args.percentile_basis]
            else:
                thr_min = pcts['front']  # no cells behind the front: fall back
            label = f'P{args.percentile:g}, {args.percentile_basis}'
        else:
            cand = np.unique(np.concatenate([damage[fr], [0.0, 1.0 + 1e-9]]))
            # loss is non-increasing in threshold; binary search the smallest feasible one
            lo, hi = 0, len(cand) - 1
            while lo < hi:
                mid = (lo + hi) // 2
                loss = area[simulate_retreat(cand[mid], dyn_float, damage, area, adj,
                                             exposed0, in_reg)].sum() / fl_area
                if loss < args.max_area_loss:
                    hi = mid
                else:
                    lo = mid + 1
            thr_min = cand[lo]
            label = f'<{100 * args.max_area_loss:g}% area loss'
        thr = min(thr_min + args.growth_margin, 1.0)
        removed = simulate_retreat(thr, dyn_float, damage, area, adj, exposed0, in_reg)
        loss = area[removed].sum() / fl_area
        p10, p50, p90 = np.percentile(damage[fr], [10, 50, 90])
        print(f'{name:45s} {fr.sum():6d} {p10:6.3f} {p50:6.3f} {p90:6.3f} '
              f'{thr:6.3f} {100 * loss:6.2f}')
        result[name] = (in_reg, thr)
        plot_region(name, geom, thr, loss, removed, in_reg, front, ice, floating, damage, x, y,
                    dc, args.fig_dir, fr_idx, near_idx, 1.e3 * args.near_dist_km, label, args.min_damage)

    out, default = combine_region_values(result, ncells, overlap='max')  # overlap -> largest
    print(f'Default (median over regions) for cells outside regions / uncalibrated: {default:.3f}')

    write_calibration_file(
        args.output, 'damageCalvingThreshold', 'damage threshold for crevasse-depth calving',
        None, out, ds,
        f'Created by calibrate_damage_threshold_for_cd_calving.py from {args.input}; '
        f'criterion={args.criterion}, max_area_loss={args.max_area_loss}, '
        f'percentile={args.percentile}, basis={args.percentile_basis}, '
        f'growth_margin={args.growth_margin}')
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
