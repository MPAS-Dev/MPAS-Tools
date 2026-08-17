#!/usr/bin/env python3
"""
Convert a MALI Weertman basal-friction initial condition to parameters
for Albany's Regularized Coulomb friction law.

Computes
--------
1. Downs & Johnson-style effective pressure N
2. Area-weighted optimal scalar C:

       C = integral[ uc**q * mu / N dA ] / integral[dA]

3. Albany bed-roughness/Lambda field:

       Lambda = uc / (A * N**n)

The output is a copy of the input MALI initial-condition file with
`Lambda` and optionally `effectivePressure` added.

Notes
-----
Albany's RC law is

    beta = C * N * |u|^(q-1) /
           (|u| + Lambda * A * N^n)^q

so Lambda * A * N^n has units of velocity.

All quantities supplied here must use a mutually consistent unit system.
For typical MALI/Albany configurations:
    velocity : m yr^-1
    N        : check whether the Albany interface expects Pa or kPa
    A        : consistent with N and yr
"""

import argparse
import shutil

import numpy as np
import xarray as xr


def downs_johnson_effective_pressure(
        thickness,
        bed,
        min_fraction_overburden,
        length_scale,
        rho_i=910.0,
        rho_w=1028.0,
        gravity=9.80616):
    """
    Reproduce Albany's Downs & Johnson-style effective pressure.

    Parameters
    ----------
    thickness : ndarray
        Ice thickness H [m].
    bed : ndarray
        Bed elevation b [m], positive above sea level.
    min_fraction_overburden : float
        Albany "Minimum Fraction Overburden Pressure".
    length_scale : float
        Albany "Length Scale Factor".

        IMPORTANT: use the same units as `bed`. If the Albany parameter
        is entered in km but bed is in m here, convert it to meters before
        calling this function.
    rho_i : float
        Ice density [kg m^-3].
    rho_w : float
        Water density [kg m^-3].
    gravity : float
        Gravitational acceleration [m s^-2].

    Returns
    -------
    N : ndarray
        Effective pressure [Pa].
    """
    H = np.asarray(thickness, dtype=np.float64)
    b = np.asarray(bed, dtype=np.float64)

    # Albany:
    # f_p = 1 / (1 + exp(-bed / length_scale))
    #
    # Use a numerically stable clipped argument.
    arg = np.clip(b / length_scale, -700.0, 700.0)
    fp = 1.0 / (1.0 + np.exp(-arg))

    overburden_term = min_fraction_overburden * rho_i * H * fp
    marine_term = (1.0 - fp) * np.maximum(-rho_w * b, 0.0)

    N = gravity * np.maximum(
        rho_i * H - (overburden_term + marine_term),
        0.0
    )

    return N


def area_weighted_optimal_C(mu, N, area, uc, q, mask):
    """
    C = integral(uc^q * mu / N dA) / integral(dA)
    """
    valid = (
        mask
        & np.isfinite(mu)
        & np.isfinite(N)
        & np.isfinite(area)
        & (N > 0.0)
        & (area > 0.0)
    )

    if not np.any(valid):
        raise ValueError("No valid grounded cells available for C calculation.")

    integrand = (uc ** q) * mu[valid] / N[valid]

    C = np.sum(area[valid] * integrand) / np.sum(area[valid])

    return C, valid


def main():
    parser = argparse.ArgumentParser(
        description="Convert MALI Weertman friction IC to Regularized Coulomb."
    )

    parser.add_argument("input", help="Input MALI initial-condition NetCDF file")
    parser.add_argument("output", help="Output NetCDF file")

    parser.add_argument(
        "--critical-velocity", "--uc",
        type=float,
        required=True,
        help="Critical velocity u_c, e.g. in m/yr"
    )
    parser.add_argument(
        "--q",
        type=float,
        default=0.2,
        help="Sliding-law exponent q (default: 0.2)"
    )

    # Downs & Johnson / Albany parameters
    parser.add_argument(
        "--min-fraction-overburden",
        type=float,
        required=True,
        help='Albany "Minimum Fraction Overburden Pressure"'
    )
    parser.add_argument(
        "--pressure-length-scale",
        type=float,
        required=True,
        help=(
            'Albany "Length Scale Factor", converted to the same '
            "length units as bedTopography (normally m)"
        )
    )

    parser.add_argument("--rho-ice", type=float, default=910.0)
    parser.add_argument("--rho-water", type=float, default=1028.0)
    parser.add_argument("--gravity", type=float, default=9.80616)

    # Needed for Lambda = uc / (A N^n)
    parser.add_argument(
        "--flow-rate",
        type=float,
        required=True,
        help=(
            "Constant Glen flow rate A, in units consistent with "
            "critical velocity and effective pressure"
        )
    )
    parser.add_argument(
        "--glen-n",
        type=float,
        default=3.0,
        help="Glen-law exponent n (default: 3)"
    )

    # MALI field names
    parser.add_argument(
        "--mu-field",
        default="muFriction",
        help="Weertman friction field (default: muFriction)"
    )
    parser.add_argument(
        "--thickness-field",
        default="thickness",
        help="Ice thickness field (default: thickness)"
    )
    parser.add_argument(
        "--bed-field",
        default="bedTopography",
        help="Bed elevation field (default: bedTopography)"
    )
    parser.add_argument(
        "--area-field",
        default="areaCell",
        help="MPAS cell area field (default: areaCell)"
    )

    parser.add_argument(
        "--lambda-field",
        default="Lambda",
        help="Name for output bed-roughness field (default: Lambda)"
    )
    parser.add_argument(
        "--effective-pressure-field",
        default="effectivePressure",
        help="Name for diagnostic N field"
    )

    parser.add_argument(
        "--time-index",
        type=int,
        default=0,
        help="Time index for Time-dependent IC fields (default: 0)"
    )

    args = parser.parse_args()

    # -------------------------------------------------------------
    # Read IC
    # -------------------------------------------------------------
    ds = xr.open_dataset(args.input)

    required = [
        args.mu_field,
        args.thickness_field,
        args.bed_field,
        args.area_field,
    ]

    missing = [name for name in required if name not in ds]
    if missing:
        raise KeyError(
            f"Input file is missing required fields: {', '.join(missing)}"
        )

    def cell_field(name):
        """Extract nCells field, dropping Time if present."""
        da = ds[name]

        if "Time" in da.dims:
            da = da.isel(Time=args.time_index)

        values = np.asarray(da.values).squeeze()

        if values.ndim != 1:
            raise ValueError(
                f"{name} must reduce to a 1-D nCells field; "
                f"got shape {values.shape}"
            )

        return values.astype(np.float64)

    mu = cell_field(args.mu_field)
    H = cell_field(args.thickness_field)
    bed = cell_field(args.bed_field)
    area = cell_field(args.area_field)

    # -------------------------------------------------------------
    # Effective pressure
    # -------------------------------------------------------------
    N = downs_johnson_effective_pressure(
        thickness=H,
        bed=bed,
        min_fraction_overburden=args.min_fraction_overburden,
        length_scale=args.pressure_length_scale,
        rho_i=args.rho_ice,
        rho_w=args.rho_water,
        gravity=args.gravity,
    )

    # Grounded-ice test used by Albany:
    #
    #     rho_i H + rho_w b > 0
    #
    # Also exclude ice-free cells.
    grounded = (
        (H > 0.0)
        & (args.rho_ice * H + args.rho_water * bed > 0.0)
    )

    # -------------------------------------------------------------
    # Optimal C
    # -------------------------------------------------------------
    C, fit_mask = area_weighted_optimal_C(
        mu=mu,
        N=N,
        area=area,
        uc=args.critical_velocity,
        q=args.q,
        mask=grounded,
    )

    # -------------------------------------------------------------
    # Lambda / Albany Bed Roughness
    #
    # uc = Lambda * A * N^n
    # -------------------------------------------------------------
    Lambda = np.zeros_like(N)

    lambda_mask = grounded & np.isfinite(N) & (N > 0.0)

    Lambda[lambda_mask] = (
        args.critical_velocity
        / (args.flow_rate * N[lambda_mask] ** args.glen_n)
    )

    # Floating/ice-free cells are deliberately zero.
    Lambda[~lambda_mask] = 0.0

    # -------------------------------------------------------------
    # Diagnostics
    # -------------------------------------------------------------
    local_C = np.full_like(N, np.nan)
    local_C[fit_mask] = (
        args.critical_velocity ** args.q
        * mu[fit_mask]
        / N[fit_mask]
    )

    print()
    print("MALI Weertman -> Regularized Coulomb conversion")
    print("------------------------------------------------")
    print(f"Input file                    : {args.input}")
    print(f"Critical velocity, uc         : {args.critical_velocity:g}")
    print(f"Power exponent, q             : {args.q:g}")
    print(f"Glen exponent, n              : {args.glen_n:g}")
    print(f"Flow rate, A                  : {args.flow_rate:.10e}")
    print(f"Cells used in C fit           : {np.count_nonzero(fit_mask)}")
    print(f"Grounded area used            : {np.sum(area[fit_mask]):.10e}")
    print()
    print(f"Optimal C                     : {C:.16e}")
    print()
    print(
        "N range on fit domain        : "
        f"{np.nanmin(N[fit_mask]):.6e} -- "
        f"{np.nanmax(N[fit_mask]):.6e}"
    )
    print(
        "Lambda range grounded        : "
        f"{np.nanmin(Lambda[lambda_mask]):.6e} -- "
        f"{np.nanmax(Lambda[lambda_mask]):.6e}"
    )
    print(
        "Local C range                : "
        f"{np.nanmin(local_C[fit_mask]):.6e} -- "
        f"{np.nanmax(local_C[fit_mask]):.6e}"
    )
    print()

    # -------------------------------------------------------------
    # Write copy of original IC.
    # Use shutil first so unrelated variables/encoding remain intact.
    # -------------------------------------------------------------
    ds.close()
    shutil.copy2(args.input, args.output)

    out = xr.open_dataset(args.output)

    ncell_dim = ds_dims = None

    # Determine nCells dimension from areaCell.
    for dim in out[args.area_field].dims:
        if dim.lower() == "ncells":
            ncell_dim = dim
            break

    if ncell_dim is None:
        # Standard MALI name is nCells; this gives a useful fallback.
        if "nCells" in out.dims:
            ncell_dim = "nCells"
        else:
            raise ValueError("Could not identify the nCells dimension.")

    out.load()
    out.close()

    # Re-open writable through xarray and rewrite.
    # For very large production files, netCDF4.Dataset can instead be
    # used to modify the copied file in-place.
    out = xr.open_dataset(args.output).load()

    out[args.lambda_field] = xr.DataArray(
        Lambda,
        dims=(ncell_dim,),
        attrs={
            "long_name": "Albany regularized-Coulomb bed roughness Lambda",
            "description": "Lambda = u_c / (A N^n)",
        },
    )

    out[args.effective_pressure_field] = xr.DataArray(
        N,
        dims=(ncell_dim,),
        attrs={
            "long_name": "Downs-Johnson effective pressure",
            "units": "Pa",
        },
    )

    # Save conversion information globally.
    out.attrs["regularizedCoulomb_C"] = float(C)
    out.attrs["regularizedCoulomb_criticalVelocity"] = (
        float(args.critical_velocity)
    )
    out.attrs["regularizedCoulomb_q"] = float(args.q)
    out.attrs["regularizedCoulomb_GlenN"] = float(args.glen_n)
    out.attrs["regularizedCoulomb_flowRate"] = float(args.flow_rate)
    out.attrs["regularizedCoulomb_minFractionOverburden"] = (
        float(args.min_fraction_overburden)
    )
    out.attrs["regularizedCoulomb_pressureLengthScale"] = (
        float(args.pressure_length_scale)
    )

    # xarray cannot safely overwrite an open source file, so use temp.
    tmp = args.output + ".tmp"
    out.to_netcdf(tmp)
    out.close()

    shutil.move(tmp, args.output)

    print(f"Wrote converted IC: {args.output}")
    print(f"Use C = {C:.16e} in the Albany RC configuration.")


if __name__ == "__main__":
    main()
