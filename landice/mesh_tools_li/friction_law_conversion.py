#!/usr/bin/env python3
"""
Convert a MALI Weertman basal-friction initial condition to parameters
for Albany's Regularized Coulomb friction law.

Computes
--------
1. Downs & Johnson-style effective pressure N
2. Area-weighted optimal scalar C:

       C = integral[ uc**qW * mu / N dA ] / integral[dA]

3. Albany bed-roughness/Lambda field:

       Lambda = uc / (A * N**n)

The output is a copy of the input MALI initial-condition file with
`Lambda` and optionally `effectivePressure` added.

Notes
-----
The input Weertman/Power-Law friction law is

    beta = mu * N * |u|^(qW-1)

with its own "Power Exponent" qW (the exponent MALI's `muFriction`
field was calibrated with). Albany's RC law is

    beta = C * N * |u|^(qR-1) /
           (|u| + Lambda * A * N^n)^qR

with its own, independent, "Power Exponent" qR. qW and qR need not
match: qW describes the input Weertman law used to derive C, while
qR is the exponent Albany will actually use for the RC law. Per
project convention, qR is always fixed at 1/3 (see RC_POWER_EXPONENT
below) and is not user-configurable.

so Lambda * A * N^n has units of velocity.

The Glen flow-rate factor A can be obtained in one of two ways,
selected via --flow-rate-type:

- "temperature" (default): computed with Albany's "Temperature Based"
  Flow Rate Type (LandIce_FlowRate_Def.hpp), i.e. a two-branch
  Arrhenius law keyed on ice temperature, evaluated using the
  basal-most (last) vertical level of the MALI `temperature` field as
  an approximation of basal temperature. This is only valid for a
  Glen's Law n of 3 (see ALBANY_FLOW_RATE_* constants).
- "constant": a single scalar value supplied via --flow-rate is used
  for all cells (Albany Flow Rate Type: Constant).

All quantities supplied here must use a mutually consistent unit system.
For typical MALI/Albany configurations:
    velocity : m yr^-1
    N        : check whether the Albany interface expects Pa or kPa
    A        : Pa^-3 s^-1 (Albany's Temperature Based flow rate is SI)
"""

import argparse
import shutil

import numpy as np
import xarray as xr

# Albany's Regularized Coulomb law always uses a fixed Power Exponent
# of 1/3. This is independent of the input Weertman/Power-Law exponent
# (--weertman-q), which is used only to derive the optimal C.
RC_POWER_EXPONENT = 1.0 / 3.0

# Seconds per year, matching Albany's own hardcoded conversion factor
# exactly (LandIce_BasalFrictionCoefficient_Def.hpp: "secsInYr = 365 *
# 24 * 3600"). Albany's Regularized Coulomb evaluator computes the
# Glen flow rate A in SI units (Pa^-3 s^-1, per second) but solves for
# velocity in m/yr; this factor converts between the two so that a
# critical velocity supplied in m/yr and a flow rate in Pa^-3 s^-1
# combine correctly when deriving Lambda (bedRoughnessRC).
SECONDS_PER_YEAR = 365.0 * 24.0 * 3600.0

# Albany's internal effective pressure representation (used e.g. by
# its "Hydrostatic" Effective Pressure Type) is computed from bed/
# thickness fields that MALI's coupling interface has already divided
# by 1000 (m -> km) before Albany ever sees them
# (Interface_velocity_solver.cpp: "unit_length = 1000"). Combined with
# SI density/gravity, this makes Albany's internal effective-pressure
# number equal to the physical pressure in Pa divided by 1000 (i.e.
# numerically "kPa", matching Albany's own documentation: "Effective
# Pressure [kPa]"). C must be derived against this same kPa-scaled N
# so that it is dimensionally consistent with Albany's own
# "beta = C * N * |u|^(q-1)" evaluation at runtime.
ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT = 1000.0

# Albany's "Temperature Based" Flow Rate Type constants, reproduced
# exactly from LandIce_FlowRate_Def.hpp. These are only valid for a
# Glen's Law n of 3 (arrmlh/arrmll are in Pa^-3 s^-1).
ALBANY_FLOW_RATE_ACTENH = 1.39e5      # [J mol-1]
ALBANY_FLOW_RATE_ACTENL = 6.0e4       # [J mol-1]
ALBANY_FLOW_RATE_GASCON = 8.314       # [J mol-1 K-1]
ALBANY_FLOW_RATE_SWITCHING_T = 263.15  # [K]
ALBANY_FLOW_RATE_ARRMLH = 1.733e3     # [Pa-3 s-1]
ALBANY_FLOW_RATE_ARRMLL = 3.613e-13   # [Pa-3 s-1]


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


def albany_temperature_based_flow_rate(temperature):
    """
    Reproduce Albany's "Temperature Based" Flow Rate Type exactly
    (LandIce_FlowRate_Def.hpp, TEMPERATURE_BASED case):

        A(T) = arrmll * exp(-actenl / (gascon * T))   if T < switchingT
        A(T) = arrmlh * exp(-actenh / (gascon * T))   otherwise

    Parameters
    ----------
    temperature : ndarray
        Ice temperature [K]. Only valid for a Glen's Law n of 3.

    Returns
    -------
    A : ndarray
        Glen flow-rate factor [Pa^-3 s^-1].
    """
    T = np.asarray(temperature, dtype=np.float64)

    A_low = ALBANY_FLOW_RATE_ARRMLL * np.exp(
        -ALBANY_FLOW_RATE_ACTENL / (ALBANY_FLOW_RATE_GASCON * T)
    )
    A_high = ALBANY_FLOW_RATE_ARRMLH * np.exp(
        -ALBANY_FLOW_RATE_ACTENH / (ALBANY_FLOW_RATE_GASCON * T)
    )

    return np.where(T < ALBANY_FLOW_RATE_SWITCHING_T, A_low, A_high)


def area_weighted_optimal_C(mu, N, area, uc, q, mask):
    """
    C = integral(uc^q * mu / N dA) / integral(dA)

    `q` here is the input Weertman/Power-Law exponent (qW), not the
    Regularized Coulomb exponent.

    IMPORTANT: `N` must already be expressed in the same units Albany
    will actually use at runtime for its own internal effective
    pressure (numerically equal to physical Pa / 1000, i.e. "kPa";
    see ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT), not raw SI Pascals.
    Passing raw-Pa N here would make C inconsistent with how Albany
    multiplies C against its own internal N at runtime.
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
        description="Convert MALI Weertman friction IC to Regularized Coulomb.",
        # Disable prefix-abbreviation matching. Without this, obsolete or
        # mistyped flags (e.g. a leftover "--flow-rate VALUE" from before
        # flow rate became Temperature Based-only) can silently be matched
        # as an unambiguous prefix of another option (e.g.
        # "--flow-rate-field"), corrupting its value instead of erroring.
        allow_abbrev=False,
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
        "--weertman-q", "--q",
        dest="weertman_q",
        type=float,
        default=0.2,
        help=(
            "Input Weertman/Power-Law sliding exponent qW, used to "
            "derive the optimal C (default: 0.2). This is independent "
            "of the Regularized Coulomb law's Power Exponent, which is "
            f"always {RC_POWER_EXPONENT:g} (see RC_POWER_EXPONENT)."
        )
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

    # Needed for Lambda = uc / (A N^n).
    parser.add_argument(
        "--flow-rate-type",
        choices=["temperature", "constant"],
        default="temperature",
        help=(
            'How to obtain the Glen flow rate A (default: temperature). '
            '"temperature" computes a per-cell A from ice temperature '
            "via Albany's Temperature Based Flow Rate Type "
            '(see --temperature-field); "constant" uses a single '
            "scalar value supplied via --flow-rate for all cells."
        )
    )
    parser.add_argument(
        "--temperature-field",
        default="temperature",
        help=(
            "MALI ice temperature field [K] (default: temperature). "
            "Used to compute the Glen flow rate A via Albany's "
            "Temperature Based Flow Rate Type when "
            "--flow-rate-type=temperature. If the field has a "
            "nVertLevels dimension, the last (basal-most) level is "
            "used as an approximation of basal temperature."
        )
    )
    parser.add_argument(
        "--flow-rate",
        type=float,
        default=None,
        help=(
            "Constant Glen flow rate A [Pa^-3 s^-1], required when "
            "--flow-rate-type=constant. Ignored otherwise."
        )
    )
    parser.add_argument(
        "--glen-n",
        type=float,
        default=3.0,
        help=(
            "Glen-law exponent n (default: 3). Albany's Temperature "
            "Based Flow Rate Type constants are only valid for n=3; "
            "a warning is issued if a different value is used with "
            "--flow-rate-type=temperature."
        )
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
        "--velocity-x-field",
        default="uReconstructX",
        help=(
            "MALI x-velocity field [m s^-1] (default: uReconstructX). "
            "Used, together with --velocity-y-field, to compute the "
            "basal sliding speed (last nVertInterfaces level) that "
            "Lambda is solved against. Units are assumed m/s, matching "
            "standard MALI output/restart files."
        )
    )
    parser.add_argument(
        "--velocity-y-field",
        default="uReconstructY",
        help="MALI y-velocity field [m s^-1] (default: uReconstructY)"
    )

    parser.add_argument(
        "--lambda-field",
        default="bedRoughnessRC",
        help=(
            "Name for output Regularized Coulomb bed-roughness field "
            "(default: bedRoughnessRC)"
        )
    )
    parser.add_argument(
        "--effective-pressure-field",
        default="effectivePressure",
        help="Name for diagnostic N field"
    )
    parser.add_argument(
        "--flow-rate-field",
        default="flowRateA",
        help="Name for diagnostic Glen flow-rate A field"
    )

    parser.add_argument(
        "--time-index",
        type=int,
        default=0,
        help="Time index for Time-dependent IC fields (default: 0)"
    )

    args = parser.parse_args()

    if args.flow_rate_type == "constant":
        if args.flow_rate is None:
            parser.error(
                "--flow-rate is required when --flow-rate-type=constant"
            )
    else:
        if args.flow_rate is not None:
            parser.error(
                "--flow-rate is only used with --flow-rate-type=constant "
                "(got --flow-rate-type=temperature)"
            )
        if args.glen_n != 3.0:
            print(
                "WARNING: Albany's Temperature Based Flow Rate Type "
                "constants (arrmlh/arrmll) are only valid for a Glen's "
                f"Law n of 3; --glen-n was set to {args.glen_n:g}. The "
                "computed flow rate A will be physically inconsistent "
                "with this exponent."
            )

    # -------------------------------------------------------------
    # Read IC
    # -------------------------------------------------------------
    ds = xr.open_dataset(args.input)

    required = [
        args.mu_field,
        args.thickness_field,
        args.bed_field,
        args.area_field,
        args.velocity_x_field,
        args.velocity_y_field,
    ]
    if args.flow_rate_type == "temperature":
        required.append(args.temperature_field)

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

    def basal_cell_field(name):
        """
        Extract nCells field from a (Time, nCells, nVertLevels),
        (Time, nCells, nVertInterfaces), (nCells, nVertLevels), or
        (nCells, nVertInterfaces) field, taking the last vertical
        level/interface as an approximation of the basal-most value.
        """
        da = ds[name]

        if "Time" in da.dims:
            da = da.isel(Time=args.time_index)

        vert_dims = [d for d in da.dims if d.lower().startswith("nvert")]
        if vert_dims:
            da = da.isel({vert_dims[0]: -1})

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
    # Basal sliding speed, from the input file's basal-most
    # (last nVertInterfaces level) horizontal velocity components.
    # This is the *actual* current sliding speed used to solve the
    # Regularized Coulomb law for Lambda below (as opposed to an
    # assumed critical velocity).
    # -------------------------------------------------------------
    uX = basal_cell_field(args.velocity_x_field)
    uY = basal_cell_field(args.velocity_y_field)
    # Convert from m/s (standard MALI units) to m/yr, matching
    # Albany's internal u_norm convention (see
    # LandIce_BasalFrictionCoefficient_Def.hpp: "Sliding Velocity
    # Regularization [m yr^-1]").
    speed = np.sqrt(uX ** 2 + uY ** 2) * SECONDS_PER_YEAR

    # -------------------------------------------------------------
    # Glen flow-rate factor A.
    # -------------------------------------------------------------
    if args.flow_rate_type == "temperature":
        basal_temperature = basal_cell_field(args.temperature_field)
        A = albany_temperature_based_flow_rate(basal_temperature)
    else:
        basal_temperature = None
        A = np.full_like(H, args.flow_rate)

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

    # Albany's own internal effective pressure (e.g. "Hydrostatic"
    # Effective Pressure Type) is computed from bed/thickness fields
    # that MALI's coupling interface has already divided by 1000 (m ->
    # km) before Albany sees them, which numerically makes Albany's
    # internal N equal to physical N[Pa] / 1000. C must be derived
    # against this same scale for consistency with Albany's own
    # "beta = C * N * |u|^(q-1)" evaluation.
    N_albany = N / ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT

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
        N=N_albany,
        area=area,
        uc=args.critical_velocity,
        q=args.weertman_q,
        mask=grounded,
    )

    # -------------------------------------------------------------
    # Lambda / Albany Bed Roughness
    #
    # Rather than assuming the sliding speed equals a prescribed
    # critical velocity, Lambda is now solved for exactly, per cell,
    # by requiring that Albany's Regularized Coulomb law reproduce
    # the *same basal shear stress* Tau_b that the original
    # Weertman/Power-Law would produce at the cell's actual current
    # sliding speed (from --velocity-x-field/--velocity-y-field).
    #
    # NOTE: MALI's Weertman sliding law (which muFriction was
    # calibrated for) has NO effective pressure term:
    #
    #   Weertman:            beta_W  = mu * u^(qW-1)
    #                        Tau_b   = beta_W * u = mu * u^qW
    #
    #   Regularized Coulomb (LandIce_BasalFrictionCoefficient_Def.hpp):
    #                        beta_RC = C * N * u^(p-1)
    #                                  / (u + Lambda*scaling*A*N^n)^p
    #                        Tau_b   = beta_RC * u
    #                                = C * N * u^p
    #                                  / (u + Lambda*scaling*A*N^n)^p
    #
    #   Setting Tau_b_RC == Tau_b_W and solving for Lambda:
    #
    #     u + Lambda*scaling*A*N^n = u * (C*N / Tau_b_W)^(1/p)
    #                              = u * (C*N / (mu * u^qW))^(1/p)
    #
    #     Lambda = u * [(C*N / (mu * u^qW))^(1/p) - 1]
    #              / (SECONDS_PER_YEAR * A * N^n)
    #
    # This is exactly the same premise used by area_weighted_optimal_C
    # above (which likewise assumes Tau_b_W = mu * uc^qW with no N
    # term, matched against the RC Coulomb limit C*N).
    #
    # where u is in m/yr, A is in Pa^-3 s^-1, N (raw Pa) is used here
    # exactly as in the previous uc-based derivation (Albany's
    # internal km/kPa/yr "scaling" factor reduces to the plain
    # SECONDS_PER_YEAR factor once Lambda is expressed in meters and N
    # in Pa). N_albany (the kPa-equivalent convention) is used for the
    # "C*N" Coulomb-limit term, matching how C was itself fit.
    #
    # Because Albany's Regularized Coulomb law can never produce a
    # shear stress above the Coulomb limit C*N (attained only in the
    # Lambda -> 0 limit), cells where the Weertman law's Tau_b at the
    # current speed already meets or exceeds C*N have no valid
    # (non-negative) solution for Lambda; these are set to 0 (maximal
    # Coulomb sliding) and reported below.
    # -------------------------------------------------------------
    Lambda = np.zeros_like(N)

    valid_speed = (
        grounded
        & np.isfinite(N) & (N > 0.0)
        & np.isfinite(mu) & (mu > 0.0)
        & np.isfinite(speed) & (speed > 0.0)
    )

    tau_b_weertman = np.full_like(N, np.nan)
    tau_b_weertman[valid_speed] = (
        mu[valid_speed] * speed[valid_speed] ** args.weertman_q
    )

    stress_ratio = np.full_like(N, np.nan)
    stress_ratio[valid_speed] = (
        (C * N_albany[valid_speed]) / tau_b_weertman[valid_speed]
    )

    lambda_mask = valid_speed & np.isfinite(stress_ratio) & (stress_ratio > 1.0)

    Lambda[lambda_mask] = (
        speed[lambda_mask]
        * (stress_ratio[lambda_mask] ** (1.0 / RC_POWER_EXPONENT) - 1.0)
        / (SECONDS_PER_YEAR * A[lambda_mask] * N[lambda_mask] ** args.glen_n)
    )

    n_unreachable = int(np.count_nonzero(valid_speed & ~lambda_mask))
    if n_unreachable > 0:
        print(
            f"WARNING: {n_unreachable} grounded cells have a Weertman "
            "basal shear stress (mu*u^qW) at the current sliding "
            "speed that meets or exceeds the Coulomb limit C*N; "
            "Lambda set to 0.0 (maximal Coulomb sliding) at these "
            "cells."
        )

    # Floating/ice-free/invalid-speed/unreachable-stress cells are
    # deliberately left at zero (see Lambda initialization above).

    # -------------------------------------------------------------
    # Diagnostics
    # -------------------------------------------------------------
    local_C = np.full_like(N, np.nan)
    local_C[fit_mask] = (
        args.critical_velocity ** args.weertman_q
        * mu[fit_mask]
        / N_albany[fit_mask]
    )

    print()
    print("MALI Weertman -> Regularized Coulomb conversion")
    print("------------------------------------------------")
    print(f"Input file                    : {args.input}")
    print(f"Critical velocity, uc         : {args.critical_velocity:g}")
    print(f"Weertman power exponent, qW   : {args.weertman_q:g}")
    print(f"RC power exponent, qR         : {RC_POWER_EXPONENT:g}")
    print(f"Glen exponent, n              : {args.glen_n:g}")
    print(f"Flow rate type                : {args.flow_rate_type}")
    print(f"Flow rate A range              : {np.min(A):.10e} -- {np.max(A):.10e}")
    if basal_temperature is not None:
        print(
            "Basal temperature range        : "
            f"{np.min(basal_temperature):.6f} -- "
            f"{np.max(basal_temperature):.6f} K"
        )
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
        "Basal sliding speed range     : "
        + (
            f"{np.nanmin(speed[valid_speed]):.6e} -- "
            f"{np.nanmax(speed[valid_speed]):.6e} m/yr"
            if np.any(valid_speed) else "n/a (no valid cells)"
        )
    )
    print(
        f"Cells with valid Lambda solve  : {int(np.count_nonzero(lambda_mask))} "
        f"/ {int(np.count_nonzero(grounded))} grounded"
    )
    print(
        "Lambda range grounded        : "
        + (
            f"{np.nanmin(Lambda[lambda_mask]):.6e} -- "
            f"{np.nanmax(Lambda[lambda_mask]):.6e}"
            if np.any(lambda_mask) else "n/a (no valid cells)"
        )
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

    # The Weertman muFriction field is not used by the Regularized
    # Coulomb law; drop it from the converted IC.
    if args.mu_field in out:
        out = out.drop_vars(args.mu_field)

    out[args.lambda_field] = xr.DataArray(
        Lambda,
        dims=(ncell_dim,),
        attrs={
            "long_name": "Albany regularized-Coulomb bed roughness Lambda",
            "description": (
                "Lambda solved exactly so that the Regularized Coulomb "
                "law reproduces the Weertman law's basal shear stress "
                "(mu*u^qW, no effective-pressure term) at the cell's "
                "actual current sliding speed u (from velocity-x/y-"
                "field, last nVertInterfaces level): Lambda = u * "
                "[(C*N/(mu*u^qW))^(1/qR) - 1] / (SECONDS_PER_YEAR * A "
                "* N^n), with u in m/yr, A in Pa^-3 s^-1, N in Pa, "
                "matching Albany's internal secsInYr scaling in "
                "LandIce_BasalFrictionCoefficient_Def.hpp"
            ),
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

    if args.flow_rate_type == "temperature":
        flow_rate_long_name = (
            "Glen flow-rate factor A from Albany's Temperature "
            "Based Flow Rate Type, evaluated at basal temperature"
        )
        albany_flow_rate_type = "Temperature Based"
    else:
        flow_rate_long_name = (
            "Glen flow-rate factor A, constant value supplied by the "
            "user (Albany Flow Rate Type: Constant)"
        )
        albany_flow_rate_type = "Constant"

    out[args.flow_rate_field] = xr.DataArray(
        A,
        dims=(ncell_dim,),
        attrs={
            "long_name": flow_rate_long_name,
            "units": "Pa-3 s-1",
        },
    )

    # Save conversion information globally.
    out.attrs["regularizedCoulomb_C"] = float(C)
    out.attrs["regularizedCoulomb_criticalVelocity"] = (
        float(args.critical_velocity)
    )
    out.attrs["regularizedCoulomb_q"] = float(RC_POWER_EXPONENT)
    out.attrs["weertman_q"] = float(args.weertman_q)
    out.attrs["regularizedCoulomb_GlenN"] = float(args.glen_n)
    out.attrs["regularizedCoulomb_flowRateType"] = albany_flow_rate_type
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

    if args.flow_rate_type == "constant":
        flow_rate_yaml_lines = (
            f"        Flow Rate Type: Constant\n"
            f"        Flow Rate: {args.flow_rate:.16e}\n"
        )
    else:
        flow_rate_yaml_lines = (
            f"        Flow Rate Type: Temperature Based\n"
        )

    print()
    print("Suggested Albany YAML section")
    print("-----------------------------")
    print(
    f"""
    LandIce BCs:
      Basal Friction Coefficient:
        Type: Regularized Coulomb
        Coulomb Friction Coefficient: {C:.16e}
        Power Exponent: {RC_POWER_EXPONENT:.16e}
{flow_rate_yaml_lines}
        Bed Roughness Type: Field
        Effective Pressure Type: Hydrostatic At Nodes
        Use Pressurized Bed Above Sea Level: true
        Minimum Fraction Overburden Pressure: {args.min_fraction_overburden:.16e}
        Length Scale Factor: {args.pressure_length_scale / 1000.0:.16e}
    """
    )
    print(
        "NOTE: Lambda was derived assuming Albany's \"Flow Rate Type\": "
        f"\"{albany_flow_rate_type}\" (in the Viscosity/Flow Rate section "
        "of the Albany YAML"
        + (
            ", applied to the basal ice temperature"
            if albany_flow_rate_type == "Temperature Based"
            else f", with \"Flow Rate\": {args.flow_rate:.16e}"
        )
        + "); ensure that setting is used at run time, or Lambda will be "
        "inconsistent with the actual A used by Albany."
    )


if __name__ == "__main__":
    main()
