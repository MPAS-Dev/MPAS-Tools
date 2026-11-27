#!/usr/bin/env python
"""
Calibrate the eigencalving parameter for MALI's eigencalving law.

MALI's eigencalving (Levermann et al., 2012) is a rate-based law: the calving
front retreats at a velocity
    calvingVelocity = min(K * max(0, eMax) * max(0, eMin), L)    [m/s]
where eMax, eMin are the principal strain rates, K is the eigencalving
parameter (``eigencalvingParameter``, units m s; the MALI default
3.14e16 m s is 1e9 m a) and L is ``config_calving_speed_limit``.

For each Antarctic ice-shelf region (the ``ice_shelf_extended_*`` features in
geometric_features, tag ``ais_ice_shelf_calving_calibration``) this script
chooses the K that keeps the present calving front where it is by balancing
calving against advection.  The cells used are the ones MALI applies rate-based
calving to: the dynamic floating calving-front cells (dynamic floating cells
that neighbor ice-free ocean, or neighbor a non-dynamic floating cell that
itself neighbors ice-free ocean) and, with --include-grounded (the default),
the dynamic grounded marine margin cells.  For simplicity the ice velocity
and the calving direction are both taken to be normal to the front, so on front
cell i

    advance speed   u_i = surfaceSpeed_i
    calving rate    c_i = min(K g_i, L),   g_i = max(0, eMax_i) * max(0, eMin_i)
    migration rate  m_i = u_i - c_i      (positive = advance)

and the target is a zero (weighted) mean migration rate over the region's
front cells, i.e. sum(w_i c_i) = sum(w_i u_i), with w_i = sqrt(areaCell_i) as a
proxy for the front length each cell represents (--weight length, default) or
w_i = 1 (--weight equal).  Without a speed limit this has the closed form
K = sum(w_i u_i) / sum(w_i g_i); with one, K is found numerically (it is larger,
because the capped cells calve less than K g_i).  Only the main eigencalving
relation is calibrated; the removal of thin ice and of stranded (iceberg) ice
that the MALI routine also performs is ignored.

Why a speed limit: the product of the two strain rates in a single cell is very
heavy-tailed, and the calving CFL timestep limit is set by the single largest
calving velocity on the front.  Without a limit, a fraction of a percent of
cells can have calving velocities orders of magnitude above the ice speed and
make the calving CFL limit far more restrictive than the advective one.  With
--calving-speed-limit auto (the default) L is --speed-limit-factor (default 2)
times the maximum dynamic ice speed in the input, which keeps the calving CFL
limit within about a factor of two of the advective one.  A number (in m/a)
sets L directly and "none" disables it.  MALI must be run with the same value
of config_calving_speed_limit (printed in m/s by this script); the MALI default
of 100 m/s is effectively no limit.

Because K can be dominated by a few cells with large strain rates, a bootstrap
5-95% range of K (resampling the region's front cells) and the effective number
of contributing cells are reported for each region.

Output: a netCDF file with ``eigencalvingParameter(Time, nCells)`` (m s),
constant within each region.  Where regions overlap the smallest value is used
(least calving).  Cells outside all regions, and regions without extensional
front cells, get the median of the calibrated values.  One figure per region
(front migration map and histograms of advance speed vs. calving rate) is
written to --fig-dir.

Requires numpy, scipy, netCDF4, matplotlib, shapely (>=2), pyproj and
geometric_features.  Example:
  calibrate_eigencalving_parameter.py -i output.nc -o eigencalvingParameter.nc
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import ListedColormap
from scipy.optimize import brentq

from calving_calibration_common import (
    add_region_arguments, combine_region_values, draw_base_map, finish_map,
    identify_front_cells, load_regions, neighbors_from_vertices, read_bed_topography,
    read_cell_area, region_mask, write_calibration_file)

SECONDS_PER_YEAR = 365.0 * 86400.0
# MALI's default eigencalving parameter is quoted as 1e9 m a
K_UNIT = 1.0e9 * SECONDS_PER_YEAR  # m s


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True,
                   help='MALI output file with cellMask, eMax, eMin, surfaceSpeed, '
                        'xCell, yCell, dcEdge and cellsOnVertex (and ideally '
                        'bedTopography and areaCell)')
    p.add_argument('-m', '--mesh', default=None,
                   help='optional MALI mesh file (same cells as input) providing '
                        'bedTopography and/or areaCell if the input file lacks them')
    p.add_argument('-o', '--output', default='eigencalvingParameter.nc')
    add_region_arguments(p)
    p.add_argument('--include-grounded', action=argparse.BooleanOptionalAction,
                   default=True,
                   help='also calibrate on dynamic grounded marine margin cells, '
                        'which MALI applies eigencalving to (default: yes)')
    p.add_argument('--weight', choices=['length', 'equal'], default='length',
                   help='weight of each front cell in the mean migration: '
                        'sqrt(areaCell) as a front-length proxy, or equal '
                        '(default length)')
    p.add_argument('--calving-speed-limit', default='auto',
                   help="calving velocity limit L (config_calving_speed_limit) assumed "
                        "in the calibration: 'auto' (--speed-limit-factor times the "
                        "maximum dynamic ice speed), a number in m/a, or 'none' "
                        "(default auto)")
    p.add_argument('--speed-limit-factor', type=float, default=2.0,
                   help="factor on the maximum ice speed for --calving-speed-limit "
                        "auto (default 2)")
    p.add_argument('--calving-thickness', type=float, default=0.0,
                   help='MALI config_calving_thickness (m).  MALI removes front cells '
                        'thinner than this regardless of K; it is not part of the '
                        'calibration, but if given the share of front cells below it is '
                        'reported for each region and flagged where large (default 0 = off)')
    p.add_argument('--max-cap-k-ratio', type=float, default=10.0,
                   help='with a calving speed limit, a region whose calibrated K exceeds '
                        'this multiple of its uncapped K is ill-conditioned (too few '
                        'extensional front cells, so the balance is reached only by '
                        'pushing K up until most of them are at the limit) and is treated '
                        'as uncalibrated (default 10)')
    p.add_argument('--k-range', type=float, nargs=2, metavar=('KMIN', 'KMAX'), default=None,
                   help='limit the calibrated K to this range, in m a (e.g. 1e7 5e8).  A '
                        'region whose balanced K is outside it is clamped, so its mean front '
                        'migration is no longer zero (reported as a percentage of the mean '
                        'advance speed; positive = advance).  Default: no limit')
    p.add_argument('--n-boot', type=int, default=1000,
                   help='bootstrap samples for the K uncertainty range (default 1000)')
    p.add_argument('--seed', type=int, default=0, help='bootstrap random seed')
    p.add_argument('--time-index', type=int, default=0)
    p.add_argument('--fig-dir', default='figures',
                   help='directory for per-region figures (default figures)')
    return p.parse_args()


def calving_rate(K, g, limit):
    """Calving velocity K*g (m/s), limited to `limit` (m/s) if not None."""
    c = K * g
    return c if limit is None else np.minimum(c, limit)


def solve_K(w, u, g, limit):
    """Eigencalving parameter giving zero weighted mean migration.

    Solves sum(w * min(K g, limit)) = sum(w * u).  Returns (K, feasible);
    feasible is False if no K can balance advection with the cap (then K is
    returned for the uncapped balance as a fallback).  K is nan if no cell has
    g > 0.
    """
    target = (w * u).sum()
    denom = (w * g).sum()
    if denom <= 0.0:
        return np.nan, False
    K0 = target / denom
    if limit is None:
        return K0, True
    if (w * g * K0 <= w * limit).all():   # no cell is capped at the uncapped K
        return K0, True
    if (w * limit)[g > 0].sum() <= target:  # even all-capped cannot balance
        return K0, False

    def f(K):
        return (w * calving_rate(K, g, limit)).sum() - target

    hi = K0
    for _ in range(80):
        hi *= 2.0
        if f(hi) > 0.0:
            break
    else:
        return K0, False
    return brentq(f, K0, hi, xtol=1e-13 * K0, rtol=1e-13), True


def bootstrap_range(w, u, g, limit, n_boot, rng):
    """5th and 95th percentiles of K when resampling the front cells."""
    n = len(w)
    ids = rng.integers(0, n, size=(n_boot, n))
    ks = []
    for ix in ids:
        k, ok = solve_K(w[ix], u[ix], g[ix], limit)
        if np.isfinite(k) and ok:
            ks.append(k)
    if not ks:
        return np.nan, np.nan
    return tuple(np.percentile(ks, [5, 95]))


def plot_region(name, geom, K, limit, in_reg, front, ice, floating, x, y, dc, fr, grd,
                w, u, g, fig_dir):
    """Two-panel figure: front migration rates (left) and histograms of the
    front advance speed and calving rate (right)."""
    to_ma = SECONDS_PER_YEAR
    cvel = calving_rate(K, g, limit)
    mig = (u - cvel) * to_ma
    capped = np.zeros(len(fr), bool) if limit is None else (K * g > limit)
    fig, (ax_map, ax_hist) = plt.subplots(1, 2, figsize=(17, 8),
                                          gridspec_kw={'width_ratios': [1.25, 1]})
    tri, idx, extent = draw_base_map(ax_map, geom, in_reg, front, ice, floating, x, y, dc)
    if tri is not None:
        ax_map.tripcolor(tri, np.zeros(len(idx)), vmin=0, vmax=1,
                         cmap=ListedColormap(['0.88']), shading='gouraud')
    vlim = max(np.percentile(np.abs(mig), 95), 1e-6)
    for sel, marker, label in ((~grd, 'o', 'floating front cell'),
                               (grd, 's', 'grounded marine margin cell')):
        if sel.any():
            ax_map.scatter(x[fr][sel], y[fr][sel], s=40, c=mig[sel], marker=marker,
                           cmap='RdBu', vmin=-vlim, vmax=vlim, edgecolors='k',
                           linewidths=0.3)
            ax_map.scatter([], [], s=40, c='0.5', marker=marker, edgecolors='k',
                           linewidths=0.3, label=label)
    if capped.any():
        ax_map.scatter(x[fr][capped], y[fr][capped], s=110, facecolors='none',
                       edgecolors='m', linewidths=1.2)
        ax_map.scatter([], [], s=110, facecolors='none', edgecolors='m', linewidths=1.2,
                       label='calving rate at the speed limit')
    sm = plt.cm.ScalarMappable(cmap='RdBu', norm=plt.Normalize(-vlim, vlim))
    fig.colorbar(sm, ax=ax_map, shrink=0.8,
                 label='front migration rate (m/a)  [blue: advance, red: retreat]')
    finish_map(ax_map, geom, extent)
    ax_map.scatter([], [], s=30, marker='s', c='0.6', label='grounded ice')
    ax_map.scatter([], [], s=30, marker='s', c='0.88', label='floating ice')
    ax_map.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), fontsize=8, ncol=2)
    ax_map.set_title('Front migration rate: advance speed - calving rate')

    adv, calv = u * to_ma, cvel * to_ma
    hi = max(np.percentile(np.concatenate([adv, calv]), 98), 1e-6)
    bins = np.linspace(0.0, hi, 41)
    ax_hist.hist(np.clip(adv, 0, hi), bins=bins, weights=w, density=True, alpha=0.55,
                 color='tab:blue', label=f'advance speed (surfaceSpeed), n={len(fr)}')
    ax_hist.hist(np.clip(calv, 0, hi), bins=bins, weights=w, density=True, alpha=0.55,
                 color='tab:red', label='calving rate min(K g, L)')
    ax_hist.axvline(np.average(adv, weights=w), color='k', ls='--',
                    label=f'weighted mean of each = {np.average(adv, weights=w):.0f} m/a')
    if limit is not None and limit * to_ma <= hi:
        ax_hist.axvline(limit * to_ma, color='m', ls=':',
                        label=f'speed limit L = {limit * to_ma:.0f} m/a')
    ax_hist.set_xlabel(f'rate (m/a; values above {hi:.0f} shown in the last bin)')
    ax_hist.set_ylabel('weighted probability density')
    ax_hist.set_title('Front cells: advance speed vs. calving rate')
    ax_hist.legend(fontsize=8)
    cap_txt = '' if limit is None else f', L = {limit * to_ma:.0f} m/a'
    fig.suptitle(f'{name}: eigencalvingParameter = {K:.4g} m s '
                 f'({K / K_UNIT:.4f} x 1e9 m a){cap_txt}', fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=130, bbox_inches='tight')
    plt.close(fig)


def main():
    args = parse_args()
    ds = netCDF4.Dataset(args.input)
    t = args.time_index
    ncells = len(ds.dimensions['nCells'])
    mask = ds['cellMask'][t, :].astype(int)
    emax = np.asarray(ds['eMax'][t, :], dtype=float)
    emin = np.asarray(ds['eMin'][t, :], dtype=float)
    speed = np.asarray(ds['surfaceSpeed'][t, :], dtype=float)
    thickness = np.asarray(ds['thickness'][t, :], dtype=float)
    bed = read_bed_topography(ds, args.mesh, t)
    x = np.asarray(ds['xCell'][:])
    y = np.asarray(ds['yCell'][:])
    dc = float(np.median(ds['dcEdge'][:]))
    area = read_cell_area(ds, args.mesh, dc, t)
    adj = neighbors_from_vertices(np.asarray(ds['cellsOnVertex'][:]), ncells)

    cls = identify_front_cells(mask, bed, adj)
    ice, floating = cls['ice'], cls['floating']
    float_front, grd_front = cls['front'], cls['grounded_margin']
    front = float_front | grd_front if args.include_grounded else float_front
    valid = np.isfinite(emax) & np.isfinite(emin) & np.isfinite(speed)
    if (front & ~valid).any():
        print(f'Ignoring {(front & ~valid).sum()} front cells with non-finite '
              'eMax/eMin/surfaceSpeed')
    front = front & valid
    print(f'{(front & float_front).sum()} floating calving-front cells'
          + (f' and {(front & grd_front).sum()} grounded marine margin cells'
             if args.include_grounded else '')
          + f' used for calibration; {cls["dyn_float"].sum()} dynamic floating cells')

    # calving speed limit L (m/s)
    dyn = ice & ((mask & 2) > 0)
    umax = np.nanmax(speed[dyn])
    opt = str(args.calving_speed_limit).lower()
    if opt == 'none':
        limit = None
    elif opt == 'auto':
        limit = args.speed_limit_factor * umax
    else:
        limit = float(opt) / SECONDS_PER_YEAR
    if limit is None:
        print('No calving speed limit assumed')
    else:
        print(f'Calving speed limit L = {limit:.4e} m/s = {limit * SECONDS_PER_YEAR:.0f} m/a '
              f'(max dynamic ice speed {umax * SECONDS_PER_YEAR:.0f} m/a); run MALI with '
              f'config_calving_speed_limit = {limit:.4e}')

    # calving-rate factor g = max(0, eMax) * max(0, eMin)  (s^-2) and front weights
    g_all = np.where(valid, np.maximum(0.0, emax) * np.maximum(0.0, emin), 0.0)
    w_all = np.sqrt(area) if args.weight == 'length' else np.ones(ncells)

    regions = load_regions(args.geometric_features, args.region_glob, args.region_tag,
                           args.crs)
    os.makedirs(args.fig_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    result = {}
    print(f"{'region':45s} {'front':>5s} {'grnd':>4s} {'g>0':>4s} {'<u>m/a':>7s} "
          f"{'K 1e9ma':>8s} {'K p5':>7s} {'K p95':>7s} {'n_eff':>6s} {'%capd':>5s}"
          + (f" {'%thin':>5s}" if args.calving_thickness > 0 else ''))
    for name, geom in regions.items():
        in_reg = region_mask(geom, x, y)
        fr = np.where(in_reg & front)[0]
        if len(fr) == 0:
            print(f'{name:45s} no front cells -> default later')
            result[name] = (in_reg, np.nan)
            continue
        w, u, g = w_all[fr], speed[fr], g_all[fr]
        grd = grd_front[fr] & ~float_front[fr]
        K, feasible = solve_K(w, u, g, limit)
        if not np.isfinite(K):
            print(f'{name:45s} {len(fr):5d} {grd.sum():4d} {0:4d} no extensional front cells '
                  '-> default later')
            result[name] = (in_reg, np.nan)
            continue
        K_uncapped = (w * u).sum() / (w * g).sum()
        if limit is not None and K > args.max_cap_k_ratio * K_uncapped:
            print(f'{name:45s} {len(fr):5d} {grd.sum():4d} {(g > 0).sum():4d} ill-conditioned: '
                  f'K would be {K / K_UNIT:.3g} x 1e9 m a ({K / K_uncapped:.0f}x the uncapped '
                  f'{K_uncapped / K_UNIT:.4f}) -> default later')
            result[name] = (in_reg, np.nan)
            continue
        K_balanced = K
        if args.k_range is not None:
            K = float(np.clip(K, args.k_range[0] * SECONDS_PER_YEAR,
                              args.k_range[1] * SECONDS_PER_YEAR))
        k5, k95 = bootstrap_range(w, u, g, limit, args.n_boot, rng)
        contrib = w * calving_rate(K, g, limit)
        n_eff = contrib.sum() ** 2 / (contrib ** 2).sum()
        pthin = 100.0 * (thickness[fr] < args.calving_thickness).mean()
        pcap = 0.0 if limit is None else 100.0 * (K * g > limit).sum() / max((g > 0).sum(), 1)
        print(f'{name:45s} {len(fr):5d} {grd.sum():4d} {(g > 0).sum():4d} '
              f'{np.average(u, weights=w) * SECONDS_PER_YEAR:7.0f} {K / K_UNIT:8.4f} '
              f'{k5 / K_UNIT:7.4f} {k95 / K_UNIT:7.4f} {n_eff:6.1f} {pcap:5.1f}'
              + (f' {pthin:5.1f}' + (' (thin-ice removal will dominate here)' if pthin > 25 else '')
                 if args.calving_thickness > 0 else '')
              + ('' if K == K_balanced else
                 f'  CLAMPED from {K_balanced / K_UNIT:.4f}: mean migration '
                 f'{100 * (w * (u - calving_rate(K, g, limit))).sum() / (w * u).sum():+.0f}% '
                 f'of mean advance')
              + ('' if feasible else '  (speed limit too low to balance advection: '
                 'using the uncapped K)'))
        result[name] = (in_reg, K)
        plot_region(name, geom, K, limit, in_reg, front, ice, floating, x, y, dc, fr,
                    grd, w, u, g, args.fig_dir)

    out, default = combine_region_values(result, ncells, overlap='min')  # overlap -> smallest
    print(f'Default (median over regions) for cells outside regions / uncalibrated: '
          f'{default / K_UNIT:.4f} x 1e9 m a')

    limit_txt = 'none' if limit is None else f'{limit:.6e} m/s'
    write_calibration_file(
        args.output, 'eigencalvingParameter',
        'eigencalving proportionality constant K2+-, calibrated for zero mean front '
        'migration', 'm s', out, ds,
        f'Created by calibrate_eigencalving_parameter.py from {args.input}; '
        f'weight={args.weight}, include_grounded={args.include_grounded}, '
        f'k_range={args.k_range} m a, '
        f'calving_speed_limit={limit_txt} (use as config_calving_speed_limit)')
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
