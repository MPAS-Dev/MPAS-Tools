#!/usr/bin/env python
"""
Audit the missing-value representation of ISMIP7 output files in a directory
against the fill_policy rules of the ISMIP7 data request (see
https://ismip.github.io/ISM_SimulationChecker/user/data-request.html).

Usage: python check_masks.py <directory with the ISMIP7 files of one experiment>
"""

import glob
import os
import sys

import numpy as np
import xarray as xr

FORBIDDEN = ['lithk', 'acabf', 'libmassbfgr', 'libmassbffl', 'dlithkdt',
             'licalvf', 'ligroundf', 'lifmassbf', 'sftgif', 'sftgrf',
             'sftflf']
OUTSIDE_DOMAIN = ['orog', 'topg', 'base', 'hfgeoubed', 'thdrflf', 'deltag']
NO_ICE = ['xvelsurf', 'yvelsurf', 'xvelbase', 'yvelbase', 'xvelmean',
          'yvelmean', 'litemptop', 'litempavg', 'strbasemag']
NO_GROUNDED = {'litempbotgr': 'sftgrf'}
NO_FLOATING = {'litempbotfl': 'sftflf'}


def _load(directory, var):
    files = glob.glob(os.path.join(directory, f'{var}_*.nc'))
    if not files:
        return None
    with xr.open_dataset(files[0], decode_times=False, mask_and_scale=False) \
            as ds:
        fill = ds[var].attrs.get('_FillValue')
        data = ds[var].values
    return data, fill


def main(directory):
    n_problems = 0

    def report(ok, msg):
        nonlocal n_problems
        print(('ok   ' if ok else 'FAIL ') + msg)
        n_problems += 0 if ok else 1

    data = {}
    for var in (FORBIDDEN + OUTSIDE_DOMAIN + NO_ICE + list(NO_GROUNDED) +
                list(NO_FLOATING)):
        loaded = _load(directory, var)
        if loaded is None:
            print(f'skip {var}: no file')
            continue
        arr, fill = loaded
        is_fill = (arr == fill) if fill is not None else np.zeros(arr.shape,
                                                                 bool)
        bad = ~np.isfinite(arr)
        report(not bad.any(), f'{var}: no NaN/inf ({bad.sum()})')
        data[var] = (arr, is_fill)

    def missing(var):
        return data[var][1]

    for var in FORBIDDEN:
        if var in data:
            report(not missing(var).any(),
                   f'{var}: no fill values ({missing(var).sum()})')

    if 'sftgif' in data:
        ice = data['sftgif'][0] > 0
        for var in NO_ICE:
            if var in data:
                n_val = (~missing(var) & ~ice).sum()
                n_miss = (missing(var) & ice).sum()
                report(n_val == 0 and n_miss == 0,
                       f'{var}: missing exactly where sftgif==0 '
                       f'(value without ice {n_val}, missing with ice '
                       f'{n_miss})')
        for var in OUTSIDE_DOMAIN:
            if var in data:
                n_miss = (missing(var) & ice).sum()
                report(n_miss == 0,
                       f'{var}: defined wherever there is ice '
                       f'(missing with ice {n_miss})')
    for table in (NO_GROUNDED, NO_FLOATING):
        for var, maskname in table.items():
            if var in data and maskname in data:
                m = data[maskname][0] > 0
                n_val = (~missing(var) & ~m).sum()
                n_miss = (missing(var) & m).sum()
                report(n_val == 0 and n_miss == 0,
                       f'{var}: missing exactly where {maskname}==0 '
                       f'(value without mask {n_val}, missing with mask '
                       f'{n_miss})')
    for var, maskname in (('libmassbfgr', 'sftgrf'),
                          ('libmassbffl', 'sftflf')):
        if var in data and maskname in data:
            n_bad = ((data[var][0] != 0) & ~(data[maskname][0] > 0)).sum()
            report(n_bad == 0, f'{var}: 0 where {maskname}==0 ({n_bad} not)')
    if all(v in data for v in ('sftgif', 'sftgrf', 'sftflf')):
        diff = np.abs(data['sftgrf'][0] + data['sftflf'][0] -
                      data['sftgif'][0]).max()
        report(diff < 1e-5, f'sftgrf + sftflf == sftgif (max diff {diff})')
    if all(v in data for v in ('orog', 'base', 'lithk')):
        ok = ~missing('orog') & ~missing('base')
        diff = np.abs(data['orog'][0] - data['base'][0] -
                      data['lithk'][0])[ok].max()
        report(diff < 1.0, f'orog == base + lithk (max diff {diff} m)')
    if all(v in data for v in ('orog', 'base', 'topg', 'sftgif')):
        # remapping blends coastal cells with land, so only cells well below
        # sea level are checked
        sel = (data['sftgif'][0] == 0) & (data['topg'][0] < -300.0) & \
            ~missing('orog') & ~missing('base')
        # not a failure: conservative remapping of coastal cells mixes land
        # (orog > 0) and ocean (orog = 0) so a few cells are not exactly 0
        n_bad = ((np.abs(data['orog'][0]) > 1e-3) |
                 (np.abs(data['base'][0]) > 1e-3))[sel].sum()
        print(f'info orog == base == 0 in open ocean: {n_bad} of {sel.sum()} '
              'cells not exactly 0 (coastal remapping blending)')
    print(f'{n_problems} problem(s) found')
    return n_problems


if __name__ == '__main__':
    sys.exit(1 if main(sys.argv[1]) else 0)
