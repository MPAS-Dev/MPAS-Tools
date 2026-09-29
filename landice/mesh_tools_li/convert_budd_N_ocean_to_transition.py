#!/usr/bin/env python3
"""
Convert a MALI Budd-law `muFriction` field, calibrated against the
"ocean-connection" effective pressure N, so that it instead reproduces
the same basal shear stress under the "transition" effective pressure
parameterization.

A Budd-type power law is `tau_b = mu * N * |u|^(q-1)`. Since neither
the current velocity `u` nor the exponent `q` differ between the two
effective-pressure choices, requiring the same basal shear stress at
each cell's current state gives a purely algebraic, per-cell
conversion:

    mu_new = mu_old * N_old / N_new

where N_old is ocean_connection_effective_pressure() and N_new is
effective_pressure4() (the "transition" parameterization), both
reproduced verbatim from friction_law_conversion.py. The conversion is
only applied where N_new > 0; elsewhere mu is left unchanged from the
input field.

The transition N field (N_new) is also written to the output file as
`effectivePressure`.
"""

import argparse
import os
import subprocess
import sys

import numpy as np
import xarray as xr

RHO_ICE = 910.0
RHO_WATER = 1028.0
GRAVITY = 9.80616


def effective_pressure4(
        thickness,
        bed,
        min_fraction_overburden,
        length_scale,
        rho_i=910.0,
        rho_w=1028.0,
        gravity=9.80616,
        h_ocean=25.0):
    """
    Effective pressure with a near-ocean region followed by a bounded
    transition to a prescribed inland fraction of overburden (the
    "transition" parameterization). Reproduced verbatim from
    friction_law_conversion.py.
    """
    H = np.asarray(thickness, dtype=np.float64)
    b = np.asarray(bed, dtype=np.float64)

    ice_term = rho_i * H
    ocean_term = np.maximum(-rho_w * b, 0.0)

    # Height above flotation in bed-elevation coordinates.
    height_above_flotation = np.maximum(
        b + (rho_i / rho_w) * H,
        0.0,
    )

    # Guard against division by zero at ice-free cells (H == 0, so
    # ice_term == 0); N will end up 0 there regardless of q, since N
    # is proportional to ice_term below.
    safe_ice_term = np.where(ice_term > 0.0, ice_term, 1.0)

    # Ocean-connected effective-pressure fraction.
    q_ocean = np.where(
        ice_term > 0.0,
        np.maximum(1.0 - ocean_term / safe_ice_term, 0.0),
        0.0,
    )

    # Prescribed inland effective-pressure fraction. Matches Albany's
    # actual "Minimum Fraction Overburden Pressure" convention: it is
    # subtracted from a retained-overburden fraction of 1, not applied
    # directly as a retained fraction, so N/Pice_inland =
    # 1 - min_fraction_overburden and floatation_fraction_inland =
    # min_fraction_overburden.
    q_inland = 1.0 - min_fraction_overburden

    # Ocean-connected value at the end of the fixed region.
    q_start = np.where(
        ice_term > 0.0,
        rho_w * h_ocean / safe_ice_term,
        0.0,
    )

    q_near_ocean = q_ocean

    if length_scale == 0.0:
        q = np.where(
            height_above_flotation <= h_ocean,
            q_near_ocean,
            q_inland,
        )
    else:
        distance_into_transition = np.maximum(
            height_above_flotation - h_ocean,
            0.0,
        )

        transition_q = (
            q_inland
            - (q_inland - q_start)
            * np.exp(
                -np.log(2.0)
                * distance_into_transition / length_scale
            )
        )

        q = np.where(
            height_above_flotation <= h_ocean,
            q_near_ocean,
            transition_q,
        )

    N = gravity * ice_term * q

    return N


def ocean_connection_effective_pressure(
        thickness,
        bed,
        rho_i=910.0,
        rho_w=1028.0,
        gravity=9.80616):
    """
    Reproduce Albany's "Hydrostatic Computed At Nodes" Effective
    Pressure Type with "Use Pressurized Bed Above Sea Level: false"
    (the "ocean-connection" parameterization). Reproduced verbatim
    from friction_law_conversion.py.
    """
    H = np.asarray(thickness, dtype=np.float64)
    b = np.asarray(bed, dtype=np.float64)

    marine_term = np.maximum(-rho_w * b, 0.0)

    N = gravity * np.maximum(rho_i * H - marine_term, 0.0)

    return N


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a MALI Budd-law muFriction field from the "
            "ocean-connection effective-pressure parameterization to "
            "the transition parameterization, preserving basal shear "
            "stress."
        ),
    )
    parser.add_argument("input", help="Input MALI NetCDF file")
    parser.add_argument("output", help="Output NetCDF file")
    parser.add_argument(
        "min_fraction_overburden", type=float,
        help="Transition parameterization: minimum floatation fraction inland."
    )
    parser.add_argument(
        "length_scale", type=float,
        help="Transition parameterization: transition length scale [m]."
    )
    parser.add_argument(
        "h_ocean", type=float,
        help=(
            "Transition parameterization: height above flotation [m] "
            "below which N is purely ocean-connected."
        )
    )
    args = parser.parse_args()

    ds = xr.open_dataset(args.input)

    def cell_field(name):
        da = ds[name]
        if "Time" in da.dims:
            da = da.isel(Time=0)
        return np.asarray(da.values).squeeze().astype(np.float64)

    thickness = cell_field("thickness")
    bed = cell_field("bedTopography")
    mu_old = cell_field("muFriction")

    N_old = ocean_connection_effective_pressure(
        thickness, bed, rho_i=RHO_ICE, rho_w=RHO_WATER, gravity=GRAVITY)
    N_new = effective_pressure4(
        thickness, bed,
        min_fraction_overburden=args.min_fraction_overburden,
        length_scale=args.length_scale,
        rho_i=RHO_ICE, rho_w=RHO_WATER, gravity=GRAVITY,
        h_ocean=args.h_ocean)

    undefined = (N_new <= 0.0) & (mu_old != 0.0)
    if np.any(undefined):
        print(
            f"warning: {int(undefined.sum())} cell(s) have N_new <= 0; "
            "muFriction left unchanged there (conversion undefined)",
            file=sys.stderr,
        )

    safe_N_new = np.where(N_new > 0.0, N_new, 1.0)
    mu_new = np.where(N_new > 0.0, mu_old * N_old / safe_N_new, mu_old)

    mu_da = ds["muFriction"]
    if "Time" in mu_da.dims:
        mu_new = np.broadcast_to(mu_new, mu_da.shape)
        N_new_out = np.broadcast_to(N_new, mu_da.shape)
    else:
        N_new_out = N_new
    ds["muFriction"] = xr.DataArray(
        mu_new,
        dims=mu_da.dims,
        attrs={
            **mu_da.attrs,
            "long_name": (
                "Budd-law friction coefficient, converted from the "
                "ocean-connection to the transition effective-pressure "
                "parameterization (see convert_budd_N_ocean_to_transition.py)"
            ),
        },
    )
    ds["effectivePressure"] = xr.DataArray(
        N_new_out,
        dims=mu_da.dims,
        attrs={
            "long_name": (
                "Effective pressure N from the transition "
                "parameterization (see convert_budd_N_ocean_to_transition.py)"
            ),
            "units": "Pa",
        },
    )
    ds.load()

    # Writing NETCDF3 directly via xarray/netCDF4 is very slow for
    # large files. Instead, write fast as NETCDF4 to a temp file, then
    # use nco's `ncks` to convert to CDF5 (64-bit data, "pnetcdf")
    # format, which MPAS requires.
    tmp_nc4 = args.output + ".nc4.tmp"
    ds.to_netcdf(tmp_nc4, format="NETCDF4")
    ds.close()

    try:
        subprocess.run(
            ["ncks", "-O", "--fl_fmt=64bit_data", tmp_nc4, args.output],
            check=True,
        )
    finally:
        os.remove(tmp_nc4)

    print(f"Wrote converted IC: {args.output}")


if __name__ == "__main__":
    main()
