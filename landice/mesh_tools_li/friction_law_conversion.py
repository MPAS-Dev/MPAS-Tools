#!/usr/bin/env python3
"""
Convert a MALI Weertman/Budd basal-friction initial condition to
parameters for Albany's Regularized Coulomb friction law.

Two methods are available, selected via --method:

- "stress-match-fit" (default): computes
  1. Downs & Johnson-style effective pressure N
  2. Area-weighted optimal scalar C:

         C = integral[ uc**qW * mu * N_source / N dA ] / integral[dA]

  3. Albany bed-roughness/Lambda field:

         Lambda = uc / (A * N**n)

  Lambda is solved per-cell so that the RC law matches the source
  law's stress outside the fast-flowing (uc) region; C is a single
  scalar written uniformly to every cell.

- "transition-velocity": given a fixed transition velocity u0
  (--transition-velocity), computes, for every grounded cell,

         Lambda = u0 / (A * N**n)
         C = tau_b * (ub + u0)**(1/3) / (N * ub**(1/3))

  (ub = current sliding speed, tau_b = source-law shear stress
  mu*N_source*ub**qW), so both Lambda and a spatially-varying C
  exactly reproduce the source law's stress at every cell's current
  speed, with no fast/slow-region split. See
  solve_transition_velocity().

In both cases N can be computed with any of the
--effective-pressure-type options below. The *input* friction law's
own effective pressure N_source (see the "Notes" section below and
--source-effective-pressure-type) is independent of N and may use a
different method/parameters entirely.

MALI extends the input mesh by one cell around its boundary when
building the FEM mesh handed to Albany, assigning a fixed minimum
thickness to that extended ring, so the ice thickness/bed elevation
Albany actually sees can differ from the plain MALI initial-condition
file at/near the domain margin. This script always uses MALI's
`config_write_albany_ascii_mesh` ascii output for thickness/bed
instead: `thickness.ascii`, `bed_topography.ascii`, and
`mpas_cellID.ascii` are required to exist in `--ascii-mesh-dir`
(default: '.') and are used for every downstream calculation
(effective pressure, grounded/floating classification,
grounding-line/terminus identification, diagnostics, and plots).

The output is a copy of the input MALI initial-condition file with
`Lambda` (stored in the bed-roughness field, `--lambda-field`,
default `bedRoughnessRC`) added, the input Weertman/Budd `muFriction`
field (`--mu-field`) overwritten with the RC coefficient C (a
spatially-uniform scalar for --method=stress-match-fit, spatially
varying for --method=transition-velocity), and optionally
`effectivePressure`/`effectivePressureSource` added.

Notes
-----
The input friction law (a generalized power law, of which classic
Weertman and Budd are special cases) is

    beta = mu * N_source^1 * |u|^(qW-1)

with its own "Power Exponent" qW (the exponent MALI's `muFriction`
field was calibrated with), and its own effective pressure
`N_source`, computed independently of Albany's own RC effective
pressure N (see --source-effective-pressure-type below). Classic
Weertman is the special case where `N_source` is held at a spatially
uniform constant of 1.0 (--source-effective-pressure-type=constant,
the default); a genuine Budd-type law is obtained by selecting
"downs-johnson", "ocean-connection", "transition", or "field"
instead, so that `N_source` varies per cell like a real effective
pressure. The exponent on `N_source` is always 1 (not
user-configurable) -- this script does not support Budd variants
with a separate pressure exponent. Albany's RC law is

    beta = C * N * |u|^(qR-1) /
           (|u| + Lambda * A * N^n)^qR

with its own, independent, "Power Exponent" qR. qW and qR need not
match: qW describes the input source law used to derive C, while
qR is the exponent Albany will actually use for the RC law. Per
project convention, qR is always fixed at 1/3 (see RC_POWER_EXPONENT
below) and is not user-configurable.

so Lambda * A * N^n has units of velocity. Lambda is solved for in
physical meters and stored in the output `bedRoughnessRC` field
as-is, in meters, matching MALI's Registry.xml declaration of
`bedRoughnessRC` (units="m"). MALI's own Albany coupling interface
(mode_forward/Interface_velocity_solver.cpp) divides the raw,
meters-valued `bedRoughnessRC` field by 1000 before ever handing it
to Albany -- the exact same way it converts `bedTopography` and
`thickness` from meters to km -- so Albany's internal "scaling"
factor (which assumes an already-km-valued field) is satisfied
automatically by MALI's coupling layer. No additional m -> km
conversion should be applied here; doing so would double-convert the
value (dividing by 1000 twice), making Lambda 1000x too small by the
time it reaches Albany.

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
    velocity  : m yr^-1
    N         : check whether the Albany interface expects Pa or kPa
    A         : Pa^-3 s^-1 (Albany's Temperature Based flow rate is SI)
    N_source  : always expressed in the same kPa-equivalent
                convention as N_albany (physical Pa / 1000), whether
                held constant (--source-effective-pressure-type=
                constant, default value 1.0) or computed from a real
                formula/field
    mu        : units depend on --source-effective-pressure-type.
                When "constant" (the default, classic Weertman:
                N_source == a spatially uniform constant, default
                1.0), mu is in kPa * (yr/m)^qW (per a correction to
                MPAS-Tools' Registry.xml; NOT Pa * (yr/m)^qW), so
                that mu * speed^qW comes out already in the same
                kPa scale as N_albany. For any other
                --source-effective-pressure-type (a genuine Budd-type
                law, N_source carrying the kPa dimension itself), mu
                is instead in plain (yr/m)^qW (no kPa dimension), so
                that mu * N_source * speed^qW still comes out in kPa
                (see ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT and
                fit_coulomb_C_fast_region())
"""

import argparse
import json
import os
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

# Albany's "Bed Roughness Field Name" evaluator (Lambda,
# bedRoughnessRC here) hardcodes a "scaling" factor in
# LandIce_BasalFrictionCoefficient_Def.hpp (scaling = secsInYr *
# pow(1000, n+1); comment: "//bedRoughness in km") that assumes the
# Lambda *value it receives* is expressed in km, not meters. However,
# Albany never reads bedRoughnessRC directly from the MALI NetCDF
# file -- MALI's own coupling interface
# (mode_forward/Interface_velocity_solver.cpp: "bedRoughnessData[index]
# = bedRoughnessRC_F[iCell] / unit_length;", unit_length = 1000)
# divides the raw field by 1000 before Albany ever sees it, exactly
# as it does for bedTopography and thickness. MALI's Registry.xml
# also declares bedRoughnessRC with units="m". So the value that must
# actually be written to the output NetCDF field is Lambda in
# physical meters, with NO additional m -> km conversion applied
# here -- that conversion is already performed by MALI's coupling
# layer. (A previous version of this script divided the solved-for
# meters-based Lambda by 1000 before storing it, under the mistaken
# assumption that Albany reads the field as-is; that produced a
# double conversion, making the stored value -- and therefore the
# Lambda Albany actually uses at runtime -- 1000x too small.)

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
    transition to a prescribed inland fraction of overburden.

    Unlike downs_johnson_effective_pressure() (Albany's own internal
    "Hydrostatic"/"Hydrostatic At Nodes" formula), this is *not*
    reproduced from an Albany-internal formula. As of this writing,
    Albany has no way to consume a precomputed N field directly for
    the Regularized Coulomb law (there is no equivalent to an
    "Effective Pressure Type: Field" option here yet), so this
    parameterization cannot currently be used to run Albany at all --
    it is included for offline evaluation/comparison purposes only,
    pending upstream Albany support.

    Parameters
    ----------
    thickness : ndarray
        Ice thickness H [m].
    bed : ndarray
        Bed elevation b [m], positive above sea level.
    min_fraction_overburden : float
        Prescribed inland *minimum floatation fraction* (Pw/Pice),
        i.e. the inland effective pressure is
        N_inland = gravity * rho_i * H * (1 - min_fraction_overburden).
        Same convention as downs_johnson_effective_pressure()'s
        parameter of the same name: Albany's "Minimum Fraction
        Overburden Pressure" is subtracted from a retained-overburden
        fraction of 1 (not applied directly as a retained fraction),
        so the inland floatation fraction approaches
        min_fraction_overburden itself (its floor), not
        1 - min_fraction_overburden. Using the same convention here
        lets the same CLI value be reused across both
        effective-pressure parameterizations.
    length_scale : float
        Distance L [m, same units as `bed`/`thickness`] over which
        half the remaining difference between the near-ocean value and
        the inland target fraction is removed, moving inland from the
        point where height above flotation first exceeds `h_ocean`.
        L == 0 disables the smooth transition (a step function is used
        instead).
    rho_i : float
        Ice density [kg m^-3].
    rho_w : float
        Water density [kg m^-3].
    gravity : float
        Gravitational acceleration [m s^-2].
    h_ocean : float
        Height above flotation [m, same units as `bed`/`thickness`]
        below which the effective pressure is assumed to be set
        purely by the ocean-connected (hydrostatic) fraction, with no
        inland transition applied (default: 25.0 m; originally
        specified as 0.025 km).

    Returns
    -------
    N : ndarray
        Effective pressure [Pa].
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
    # actual "Minimum Fraction Overburden Pressure" convention (see
    # downs_johnson_effective_pressure()): it is subtracted from a
    # retained-overburden fraction of 1, not applied directly as a
    # retained fraction, so N/Pice_inland = 1 - min_fraction_overburden
    # and floatation_fraction_inland = min_fraction_overburden.
    q_inland = 1.0 - min_fraction_overburden

    # Ocean-connected value at the end of the fixed region.
    q_start = np.where(
        ice_term > 0.0,
        rho_w * h_ocean / safe_ice_term,
        0.0,
    )

    # This cap guarantees that the transition begins at or below the
    # prescribed inland value. Use np.minimum (not the Python builtin
    # min()), since q_start/q_inland are arrays/broadcastable, not
    # plain scalars.
    q_start = np.minimum(q_start, q_inland)
    # Keep the near-grounding-line branch bounded as well.
    q_near_ocean = np.minimum(q_ocean, q_inland)

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

    # Roundoff safeguard.
    q = np.clip(q, 0.0, q_inland)

    N = gravity * ice_term * q

    return N

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
        Albany "Minimum Fraction Overburden Pressure". Despite the
        name, this is *subtracted* from a retained-overburden
        fraction of 1 in Albany's actual formula (see
        LandIce_BasalFrictionCoefficient_Def.hpp), not applied
        directly as a retained fraction: far inland (bed elevation
        well above sea level), N/Pice -> 1 - min_fraction_overburden,
        i.e. the floatation fraction Pw/Pice -> min_fraction_
        overburden itself (its floor/minimum value, attained inland;
        near the grounding line the floatation fraction rises toward
        1 regardless of this parameter).
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


def ocean_connection_effective_pressure(
        thickness,
        bed,
        rho_i=910.0,
        rho_w=1028.0,
        gravity=9.80616):
    """
    Reproduce Albany's "Hydrostatic Computed At Nodes" Effective
    Pressure Type with "Use Pressurized Bed Above Sea Level: false"
    (LandIce_BasalFrictionCoefficient_Def.hpp).

    With `use_pressurized_bed` false, Albany's sigmoid bed-pressure
    fraction f_p is identically 0 everywhere (rather than a smooth
    function of bed elevation), which removes any dependence on
    "Minimum Fraction Overburden Pressure"/"Length Scale Factor" (they
    are not even read by Albany in this case) and collapses
    downs_johnson_effective_pressure()'s formula to a simple
    overburden-minus-ocean-pressure relation:

        N = max(rho_i * g * H - max(-rho_w * g * bed, 0), 0)

    i.e. full ocean (hydrostatic) water pressure wherever the bed is
    below sea level, with no inland tapering/floor -- hence
    "ocean-connection" here: N depends only on whether/how deeply the
    bed sits below sea level, not on distance from the ocean.

    Parameters
    ----------
    thickness : ndarray
        Ice thickness H [m].
    bed : ndarray
        Bed elevation b [m], positive above sea level.
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

    marine_term = np.maximum(-rho_w * b, 0.0)

    N = gravity * np.maximum(rho_i * H - marine_term, 0.0)

    return N


def compute_named_effective_pressure(
        kind,
        thickness,
        bed,
        rho_i=910.0,
        rho_w=1028.0,
        gravity=9.80616,
        min_fraction_overburden=None,
        pressure_length_scale=None,
        transition_h_ocean=None):
    """
    Dispatch to one of the three shared effective-pressure
    parameterizations by name -- "downs-johnson"
    (downs_johnson_effective_pressure()), "ocean-connection"
    (ocean_connection_effective_pressure()), or "transition"
    (effective_pressure4()) -- forwarding only the parameters each one
    actually uses. Shared by both the Regularized Coulomb law's own N
    (--effective-pressure-type) and the input source friction law's
    own N_source (--source-effective-pressure-type), so the two can
    independently select among the same three formulas (with
    independent parameters) without duplicating the dispatch logic.

    Returns
    -------
    N : ndarray
        Effective pressure [Pa].
    """
    if kind == "downs-johnson":
        return downs_johnson_effective_pressure(
            thickness=thickness,
            bed=bed,
            min_fraction_overburden=min_fraction_overburden,
            length_scale=pressure_length_scale,
            rho_i=rho_i,
            rho_w=rho_w,
            gravity=gravity,
        )
    elif kind == "ocean-connection":
        return ocean_connection_effective_pressure(
            thickness=thickness,
            bed=bed,
            rho_i=rho_i,
            rho_w=rho_w,
            gravity=gravity,
        )
    elif kind == "transition":
        return effective_pressure4(
            thickness=thickness,
            bed=bed,
            min_fraction_overburden=min_fraction_overburden,
            length_scale=pressure_length_scale,
            rho_i=rho_i,
            rho_w=rho_w,
            gravity=gravity,
            h_ocean=transition_h_ocean,
        )
    else:
        raise ValueError(f"Unsupported effective-pressure kind: {kind!r}")


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


def fit_coulomb_C_fast_region(tau_b_source, N, area, mask):
    """
    Fit a single scalar Regularized Coulomb "Coulomb Friction
    Coefficient" C by assuming that, in the fast-flowing region
    identified by `mask` (typically grounded cells with a current
    sliding speed above some critical velocity), the ice is already in
    the fully-plastic Coulomb regime of the RC law, i.e.

        Tau_b_RC = C * N

    C is chosen to be the area-weighted mean, over the fast-flowing
    region, of the per-cell "local C" implied by that relation against
    the input source law's basal shear stress at the cell's actual
    current sliding speed, `tau_b_source` (mu * N_source * speed^qW
    for a generalized power law; N_source == 1 everywhere reduces
    this to classic Weertman -- see --source-effective-pressure-type
    in main()):

        local_C = Tau_b_source / N

        C = sum(area * local_C) / sum(area)

    A straight area-weighted mean is used (rather than an N^2-weighted
    least-squares fit through the origin) so that a small number of
    high-N, high-leverage cells cannot dominate the result.

    IMPORTANT: both `tau_b_source` and `N` must already be expressed
    in the same units Albany will actually use at runtime for its own
    internal effective pressure (numerically equal to physical Pa /
    1000, i.e. "kPa"; see ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT), not
    raw SI Pascals. Passing raw-Pa N here would make C inconsistent
    with how Albany multiplies C against its own internal N at
    runtime.
    """
    valid = (
        mask
        & np.isfinite(tau_b_source)
        & np.isfinite(N)
        & np.isfinite(area)
        & (N > 0.0)
        & (area > 0.0)
        & (tau_b_source > 0.0)
    )

    if not np.any(valid):
        raise ValueError(
            "No valid fast-flowing cells (speed > critical velocity) "
            "available for C calculation."
        )

    local_C = tau_b_source[valid] / N[valid]

    C = np.sum(area[valid] * local_C) / np.sum(area[valid])

    return C, valid


def solve_transition_velocity(
    tau_b_source, speed, N, N_albany, A, glen_n, transition_velocity,
    lambda_reference_value, mu_reference_value, valid,
):
    """
    Alternative to fit_coulomb_C_fast_region() plus the per-cell
    Lambda solve: rather than fitting a single scalar C over a
    fast-flowing region and solving for Lambda elsewhere, pick a
    fixed transition velocity u0 and compute both Lambda and a
    spatially-varying Regularized Coulomb coefficient C in closed
    form, for every valid cell, such that the RC law exactly
    reproduces the input source law's basal shear stress
    (tau_b_source = mu * N_source * speed^qW; N_source == 1
    everywhere reduces this to classic Weertman -- see
    --source-effective-pressure-type in main()) at that cell's
    current sliding speed:

        Lambda = u0 / (SECONDS_PER_YEAR * A * N^n)
        C      = tau_b_source * (speed + u0)^qR
                 / (N_albany * speed^qR)

    (qR = RC_POWER_EXPONENT). This follows from substituting
    Lambda*SECONDS_PER_YEAR*A*N^n = u0 into the RC law's
    Tau_b_RC == Tau_b_source equation used elsewhere in this
    script (see the module-level Lambda-solve comment in main()) --
    i.e. u0 plays the same role as the per-cell solve's "stress
    ratio" term, but is fixed instead of solved from a separate
    scalar C fit. Unlike fit_coulomb_C_fast_region()/the per-cell
    Lambda solve, there is no fast/slow-region split and no
    ill-defined-solve case: every cell in `valid` gets an exact
    match.

    `N` (raw Pa) is used for Lambda for unit consistency with A
    (Pa^-3 s^-1) and SECONDS_PER_YEAR, exactly as in the per-cell
    Lambda solve elsewhere in this script. `N_albany` (physical N /
    ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT) is used for C, matching
    how C is fit in fit_coulomb_C_fast_region() (dimensionally
    consistent with tau_b_source).

    Cells outside `valid` (non-grounded, ice-free, or zero/invalid
    current sliding speed/N/mu/N_source) get `lambda_reference_value`
    and `mu_reference_value` respectively, since the closed-form
    expressions above are undefined there (division by zero speed,
    or an ill-defined tau_b_source/N).

    Returns
    -------
    Lambda, C : 1-D numpy arrays, same shape as `N`
    """
    Lambda = np.full_like(N, lambda_reference_value)
    Lambda[valid] = transition_velocity / (
        SECONDS_PER_YEAR * A[valid] * N[valid] ** glen_n
    )

    C = np.full_like(N, mu_reference_value)
    C[valid] = (
        tau_b_source[valid]
        * (speed[valid] + transition_velocity) ** RC_POWER_EXPONENT
        / (N_albany[valid] * speed[valid] ** RC_POWER_EXPONENT)
    )

    return Lambda, C


def cell_has_neighbor_where(test_mask, cells_on_cell, n_edges_on_cell):
    """
    For every cell, return True if at least one of its mesh neighbors
    (per MPAS `cellsOnCell`/`nEdgesOnCell` connectivity) satisfies
    `test_mask`.

    Parameters
    ----------
    test_mask : 1-D bool array (nCells,)
        Per-cell predicate to test neighbors against.
    cells_on_cell : 2-D int array (nCells, maxEdges)
        MPAS `cellsOnCell` field: 1-based neighbor cell indices, with
        0 used to pad cells with fewer than `maxEdges` edges/
        neighbors (e.g. mesh-boundary cells).
    n_edges_on_cell : 1-D int array (nCells,)
        MPAS `nEdgesOnCell` field: number of valid (non-padding)
        entries in each row of `cells_on_cell`.

    Returns
    -------
    1-D bool array (nCells,)
    """
    max_edges = cells_on_cell.shape[1]
    edge_index = np.arange(max_edges)[np.newaxis, :]
    # A neighbor slot is real only if it is within nEdgesOnCell for
    # that row *and* not the 0 padding value cellsOnCell itself uses
    # for missing neighbors (mesh-boundary cells).
    valid_slot = (
        (edge_index < n_edges_on_cell[:, np.newaxis])
        & (cells_on_cell > 0)
    )
    neighbor_index = np.where(valid_slot, cells_on_cell - 1, 0)
    neighbor_test = test_mask[neighbor_index] & valid_slot
    return np.any(neighbor_test, axis=1)


def _read_albany_ascii_field(path):
    """
    Read one of MALI's `config_write_albany_ascii_mesh` ascii fields
    (e.g. `thickness.ascii`, `bed_topography.ascii`,
    `mpas_cellID.ascii`): a first line giving the number of rows,
    followed by that many one-value-per-line rows.

    Returns
    -------
    1-D numpy array (float64), length equal to the declared row count.
    """
    values = np.loadtxt(path, dtype=np.float64)
    n_declared = int(round(values[0]))
    values = values[1:]
    if values.shape[0] != n_declared:
        raise ValueError(
            f"{path}: header declares {n_declared} rows but file has "
            f"{values.shape[0]} data rows"
        )
    return values


def load_albany_ascii_geometry(ascii_dir, n_cells):
    """
    Load the ice thickness and bed elevation actually seen by Albany
    from MALI's `config_write_albany_ascii_mesh` ascii output
    (`thickness.ascii`, `bed_topography.ascii`, `mpas_cellID.ascii`
    in `ascii_dir`), mapped onto the full MPAS `nCells` array.

    MALI extends the input mesh by one cell around its boundary when
    building the FEM mesh handed to (standalone) Albany, and may
    assign a fixed minimum ("fixed margin") thickness to that
    extended ring, so `thickness`/`bedTopography` as seen by Albany
    can differ from the plain MALI initial-condition file at/near the
    domain margin. `mpas_cellID.ascii` gives, for each row of
    `thickness.ascii`/`bed_topography.ascii`, the 1-based MPAS cell
    index that row corresponds to (see
    https://github.com/MPAS-Dev/compass/blob/main/compass/landice/
    tests/ensemble_generator/ensemble_member.py#L373-L385 for the
    same mapping convention). Values in these ascii files are in
    Albany's internal km convention and are converted to meters here
    (x1000) to match MALI's netCDF `thickness`/`bedTopography` units.

    Parameters
    ----------
    ascii_dir : str
        Directory containing `thickness.ascii`, `bed_topography.ascii`,
        and `mpas_cellID.ascii`.
    n_cells : int
        Number of cells in the MPAS mesh (`nCells`), used to validate
        `mpas_cellID.ascii` indices and to size the returned mask.

    Returns
    -------
    cell_index : 1-D int array
        0-based MPAS cell indices covered by the ascii files.
    thickness_m, bed_m : 1-D float64 arrays, same shape as
        `cell_index`
        Ice thickness and bed elevation [m], in MPAS cell order
        matching `cell_index`.
    """
    thickness_path = os.path.join(ascii_dir, "thickness.ascii")
    bed_path = os.path.join(ascii_dir, "bed_topography.ascii")
    cell_id_path = os.path.join(ascii_dir, "mpas_cellID.ascii")

    missing = [
        path for path in (thickness_path, bed_path, cell_id_path)
        if not os.path.isfile(path)
    ]
    if missing:
        raise FileNotFoundError(
            "Missing required Albany ascii mesh file(s): "
            f"{', '.join(missing)} (expected in --ascii-mesh-dir "
            f"{ascii_dir!r})"
        )

    thickness_km = _read_albany_ascii_field(thickness_path)
    bed_km = _read_albany_ascii_field(bed_path)
    cell_id = _read_albany_ascii_field(cell_id_path).astype(np.int64)

    if not (thickness_km.shape == bed_km.shape == cell_id.shape):
        raise ValueError(
            "thickness.ascii, bed_topography.ascii, and "
            "mpas_cellID.ascii must all have the same number of rows "
            f"(got {thickness_km.shape[0]}, {bed_km.shape[0]}, "
            f"{cell_id.shape[0]})"
        )

    if np.unique(cell_id).shape[0] != cell_id.shape[0]:
        raise ValueError(
            f"{cell_id_path}: cell indices are not unique"
        )

    if cell_id.min() < 1 or cell_id.max() > n_cells:
        raise ValueError(
            f"{cell_id_path}: cell indices must be in [1, {n_cells}] "
            f"(got range [{cell_id.min()}, {cell_id.max()}])"
        )

    cell_index = cell_id - 1
    return cell_index, thickness_km * 1000.0, bed_km * 1000.0


def creep_fill_extrapolate(
        values, keep_mask, fill_mask, x_cell, y_cell,
        cells_on_cell, n_edges_on_cell, method="idw",
        max_iterations=None,
):
    """
    Extrapolate `values` into every cell where `fill_mask` is True by
    repeatedly propagating values inward from `keep_mask` cells
    across MPAS mesh connectivity (`cells_on_cell`/`n_edges_on_cell`),
    a "creep fill" adapted from MPAS-Tools'
    conversion_exodus_init_to_mpasli_mesh.py (its beta/muFriction/
    stiffnessFactor extrapolation loop).

    Each iteration, every not-yet-filled `fill_mask` cell adjacent to
    at least one already-valid cell is assigned a new value derived
    from its valid neighbors only (inverse-distance weighted average,
    method="idw", or the minimum, method="min"), using the *previous*
    iteration's valid set as the source (so a single pass never
    chains through cells filled earlier in that same pass); it then
    becomes valid itself for the next iteration. This repeats until
    every `fill_mask` cell has been filled, `max_iterations` passes
    have been made, or a pass fills no new cells (a stall, meaning
    some `fill_mask` cells have no path to a `keep_mask` cell through
    other `fill_mask` cells) -- in either of the latter two cases, a
    warning is printed and any still-unfilled cells are left with
    their original `values`.

    Parameters
    ----------
    values : 1-D float array (nCells,)
        Field to extrapolate. Not modified in place; the filled
        array is returned separately.
    keep_mask : 1-D bool array (nCells,)
        True at cells whose current `values` are already valid and
        may be used as an extrapolation source.
    fill_mask : 1-D bool array (nCells,)
        True at cells whose current `values` should be discarded and
        instead derived by creep-fill extrapolation. Must be
        disjoint from `keep_mask`. Cells that are neither
        `keep_mask` nor `fill_mask` (e.g. non-grounded cells) are
        never used as a source and are left unchanged.
    x_cell, y_cell : 1-D float arrays (nCells,)
        MPAS cell-center coordinates, used for the "idw" method.
    cells_on_cell, n_edges_on_cell : see cell_has_neighbor_where()
    method : {"idw", "min"}
        "idw": inverse-distance-weighted average of valid neighbors.
        "min": minimum value among valid neighbors (matches
        MPAS-Tools' conversion_exodus_init_to_mpasli_mesh.py "min"
        extrapolation option).
    max_iterations : int, optional
        Maximum number of creep-fill passes. Default (None): no
        limit other than a stall (see above).

    Returns
    -------
    1-D float array (nCells,), same shape as `values`
    """
    if method not in ("idw", "min"):
        raise ValueError(f"Unknown creep-fill method: {method!r}")

    out = np.array(values, dtype=np.float64, copy=True)
    valid_mask = np.copy(keep_mask)
    remaining = np.copy(fill_mask)

    iteration = 0
    while np.any(remaining):
        if max_iterations is not None and iteration >= max_iterations:
            print(
                f"WARNING: creep-fill extrapolation stopped after "
                f"{iteration} iterations with "
                f"{int(np.count_nonzero(remaining))} cell(s) still "
                "unfilled; leaving their original values unchanged."
            )
            break

        newly_filled = np.zeros(out.shape, dtype=bool)

        for i_cell in np.where(remaining)[0]:
            n_edges = n_edges_on_cell[i_cell]
            neighbor_idx = cells_on_cell[i_cell, :n_edges] - 1
            neighbor_idx = neighbor_idx[neighbor_idx >= 0]

            source_idx = neighbor_idx[valid_mask[neighbor_idx]]
            if source_idx.size == 0:
                continue

            if method == "idw":
                ds = np.sqrt(
                    (x_cell[i_cell] - x_cell[source_idx]) ** 2
                    + (y_cell[i_cell] - y_cell[source_idx]) ** 2
                )
                if np.any(ds == 0.0):
                    # Degenerate (coincident) cell centers: fall back
                    # to a plain average rather than dividing by
                    # zero.
                    out[i_cell] = np.mean(out[source_idx])
                else:
                    weights = 1.0 / ds
                    out[i_cell] = (
                        np.sum(weights * out[source_idx])
                        / np.sum(weights)
                    )
            else:  # method == "min"
                out[i_cell] = np.min(out[source_idx])

            newly_filled[i_cell] = True

        if not np.any(newly_filled):
            print(
                f"WARNING: creep-fill extrapolation stalled with "
                f"{int(np.count_nonzero(remaining))} cell(s) still "
                "unfilled (no remaining cell has a valid neighbor); "
                "leaving their original values unchanged."
            )
            break

        valid_mask[newly_filled] = True
        remaining[newly_filled] = False
        iteration += 1

    return out


def load_transect(name, transects_dir):
    """
    Load a flowline transect's lon/lat coordinates directly from a
    geometric_features-style geojson file, without depending on the
    geometric_features python package (which does not yet expose
    these newer flowline transects).

    Expects `<transects_dir>/<name>/transect.geojson`, containing a
    single Feature with a LineString geometry of [lon, lat] pairs in
    degrees (geometric_features convention).

    Returns
    -------
    lon, lat : 1-D numpy arrays, in degrees.
    """
    path = os.path.join(transects_dir, name, "transect.geojson")

    with open(path) as f:
        geojson = json.load(f)

    feature = geojson["features"][0]
    geom = feature["geometry"]

    if geom["type"] != "LineString":
        raise ValueError(
            f"Transect {name!r} ({path}) has unsupported geometry "
            f"type {geom['type']!r}; only LineString is supported."
        )

    coords = np.asarray(geom["coordinates"], dtype=np.float64)
    lon = coords[:, 0]
    lat = coords[:, 1]

    return lon, lat


def project_transect(lon, lat, transformer):
    """
    Project a transect's lon/lat coordinates (degrees) into the
    planar x/y coordinate system used by the MALI mesh (via the
    supplied pyproj Transformer, e.g. EPSG:4326 -> EPSG:3031 for
    Antarctica), and compute the cumulative along-transect distance
    from the first point.

    Returns
    -------
    x, y : 1-D numpy arrays, meters, in the MALI mesh's planar CRS.
    distance : 1-D numpy array, meters, cumulative arc length along
        the transect starting from 0 at the first point.
    """
    x, y = transformer.transform(lon, lat)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    segment_length = np.hypot(np.diff(x), np.diff(y))
    distance = np.concatenate(([0.0], np.cumsum(segment_length)))

    return x, y, distance


def plot_transects(
    transect_names,
    transects_dir,
    plot_dir,
    x_cell,
    y_cell,
    fields,
    thickness,
    bed,
    rho_i,
    rho_w,
    min_fraction_overburden=None,
    floatation_fraction_label="floatation fraction",
    fit_mask=None,
    implied_c_label="implied C",
    fitted_c=None,
):
    """
    For each named transect, sample the given cell-centered fields
    (nearest-neighbor, via a KD-tree on MALI cell centers) along the
    transect and save a stacked-panel PNG plot vs. along-transect
    distance. The first panel always shows bed and ice-surface
    elevation (computed from `thickness`/`bed`); subsequent panels
    show the caller-supplied `fields`. Background shading indicates
    floating ice, ice-free ocean, and ice-free land (grounded ice is
    left unshaded). Vertical dotted gray lines mark along-transect
    distances where bed elevation crosses zero (interpolated between
    sampled points).

    Parameters
    ----------
    transect_names : list of str
        Names of subdirectories under `transects_dir`, each expected
        to contain a `transect.geojson` (geometric_features
        convention; see load_transect()).
    transects_dir : str
        Path to a geometric_features `landice/transect` directory
        (e.g. `.../geometric_features/geometric_data/landice/
        transect`).
    plot_dir : str
        Directory to write output PNGs to (created if needed).
    x_cell, y_cell : 1-D numpy arrays
        MALI mesh cell-center coordinates, meters, in the same planar
        CRS the transects will be projected into (Antarctic MALI
        meshes: EPSG:3031 polar stereographic).
    fields : dict of str -> (1-D numpy array, str, str[, (float, float)])
        Mapping of field label -> (values on MALI cells, units
        string, long_name string, optional (ymin, ymax) bound) to
        sample and plot along each transect, e.g. {"N": (N, "Pa",
        "Effective pressure")}. The optional 4th tuple element, if
        given, clamps that panel's autoscaled y-axis range to be no
        wider than (ymin, ymax) -- e.g. (0, 1) for a fraction field --
        without forcing the full range when the real data only span a
        narrower band within those bounds (which would otherwise
        squash genuinely narrow-but-valid data, e.g. 0.95-1.0, into an
        invisible sliver).
    thickness, bed : 1-D numpy arrays
        MALI cell-centered ice thickness [m] and bed elevation [m],
        used to classify each cell as grounded ice, floating ice,
        ice-free ocean, or ice-free land for background shading
        (same grounded-ice test used elsewhere in this script:
        rho_i * H + rho_w * bed > 0).
    rho_i, rho_w : float
        Ice/water density [kg m^-3], for the grounded-ice test above.
    min_fraction_overburden : float, optional
        If provided, and `fields` contains an entry keyed
        `floatation_fraction_label`, draw a horizontal reference line
        at min_fraction_overburden on that panel: the minimum
        floatation fraction (Pw/Pice) approached far inland (see
        downs_johnson_effective_pressure()/effective_pressure4()).
    floatation_fraction_label : str
        Key in `fields` identifying the floatation-fraction panel
        (default: "floatation fraction"), used to place the
        `min_fraction_overburden` reference line above.
    fit_mask : 1-D bool numpy array, optional
        MALI cell-centered mask (same shape as `thickness`/`bed`)
        marking the region used to fit the single global scalar
        Coulomb Friction Coefficient C (i.e. the "fast-flowing"/
        full-Coulomb-regime region). If provided, contiguous
        along-transect runs where this mask is True are shaded (on
        every panel) to indicate where the fit region is crossed.
    implied_c_label : str
        Key in `fields` identifying the implied-local-C panel
        (default: "implied C"). On that panel, the line is only
        drawn where `fit_mask` is True (elsewhere set to NaN so the
        line breaks), and, if `fitted_c` is given, a horizontal
        reference line at the fitted scalar C is overlaid but
        restricted to the same along-transect fit-region runs (so it
        is directly comparable to the local values plotted alongside
        it).
    fitted_c : float, optional
        The single global scalar Coulomb Friction Coefficient C, used
        to draw the reference line on the `implied_c_label` panel
        described above.
    """
    try:
        import pyproj
    except ImportError as e:
        raise ImportError(
            "Plotting transects requires the 'pyproj' package "
            "(projects transect lon/lat into the MALI mesh's planar "
            "CRS). Install it (e.g. `conda install pyproj`) or "
            "disable transect plotting with --no-plot-transects."
        ) from e

    try:
        from scipy.spatial import cKDTree
    except ImportError as e:
        raise ImportError(
            "Plotting transects requires the 'scipy' package "
            "(nearest-neighbor sampling of MALI cell fields onto "
            "transect points). Install it or disable transect "
            "plotting with --no-plot-transects."
        ) from e

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    os.makedirs(plot_dir, exist_ok=True)

    # Antarctic MALI meshes use EPSG:3031 polar stereographic (see
    # e.g. MPAS-Tools' ismip7_postprocessing/grid_and_mapping.py);
    # transects are stored in geographic lon/lat (EPSG:4326).
    transformer = pyproj.Transformer.from_crs(
        "epsg:4326", "epsg:3031", always_xy=True
    )

    tree = cKDTree(np.column_stack((x_cell, y_cell)))

    # Classify each MALI cell for background shading: grounded ice is
    # left unshaded; floating ice, ice-free ocean, and ice-free land
    # are shaded. Uses the same grounded-ice test as the rest of this
    # script (rho_i * H + rho_w * bed > 0).
    ice_free = thickness <= 0.0
    grounded = (~ice_free) & (rho_i * thickness + rho_w * bed > 0.0)
    floating = (~ice_free) & (~grounded)
    ice_free_ocean = ice_free & (bed < 0.0)
    ice_free_land = ice_free & (bed >= 0.0)

    # 0 = grounded ice (unshaded), 1 = floating ice, 2 = ice-free
    # ocean, 3 = ice-free land.
    region_class = np.zeros(thickness.shape, dtype=np.int8)
    region_class[floating] = 1
    region_class[ice_free_ocean] = 2
    region_class[ice_free_land] = 3

    shading = {
        1: ("Floating ice", "tab:blue"),
        2: ("Ice-free ocean", "tab:cyan"),
        3: ("Ice-free land", "tab:brown"),
    }

    # Upper-surface elevation, consistent with the same grounded/
    # floating classification used for shading above: grounded ice
    # surface = bed + thickness; floating ice surface = thickness *
    # (1 - rho_i/rho_w) (flotation criterion); ice-free cells show
    # bare bed/sea-level (max(bed, 0)) since there is no ice surface.
    surface = np.empty_like(thickness, dtype=np.float64)
    surface[grounded] = bed[grounded] + thickness[grounded]
    surface[floating] = thickness[floating] * (1.0 - rho_i / rho_w)
    surface[ice_free] = np.maximum(bed[ice_free], 0.0)

    # Elevation panel is always shown first, ahead of the
    # caller-supplied `fields`.
    elevation_fields = {
        "elevation": (
            {"surface": surface, "bed": bed},
            "m",
            "Bed and surface elevation",
        ),
    }
    all_fields = {**elevation_fields, **fields}

    for name in transect_names:
        lon, lat = load_transect(name, transects_dir)
        x, y, distance = project_transect(lon, lat, transformer)

        _, cell_indices = tree.query(np.column_stack((x, y)))
        distance_km = distance / 1000.0
        class_along_transect = region_class[cell_indices]
        bed_along_transect = bed[cell_indices]

        # Distances (interpolated) where bed elevation crosses zero,
        # marked with a vertical dotted line on every panel (a rough
        # proxy for the coastline/grounding-line vicinity).
        bed_sign_changes = np.flatnonzero(
            np.diff(np.sign(bed_along_transect)) != 0
        )
        zero_crossing_distances_km = []
        for i in bed_sign_changes:
            b0, b1 = bed_along_transect[i], bed_along_transect[i + 1]
            if b0 == b1:
                continue
            frac = -b0 / (b1 - b0)
            zero_crossing_distances_km.append(
                distance_km[i] + frac * (distance_km[i + 1] - distance_km[i])
            )

        fig, axes = plt.subplots(
            len(all_fields), 1, sharex=True,
            figsize=(8, 2.5 * len(all_fields))
        )
        if len(all_fields) == 1:
            axes = [axes]

        # Shade contiguous along-transect runs of each non-grounded
        # class, on every panel.
        change_indices = np.flatnonzero(
            np.diff(class_along_transect) != 0
        ) + 1
        run_starts = np.concatenate(([0], change_indices))
        run_ends = np.concatenate(
            (change_indices, [len(class_along_transect) - 1])
        )
        classes_used = set()
        for run_start, run_end in zip(run_starts, run_ends):
            cls = class_along_transect[run_start]
            if cls == 0:
                continue
            classes_used.add(cls)
            _, color = shading[cls]
            x0 = distance_km[run_start]
            x1 = distance_km[run_end]
            for ax in axes:
                ax.axvspan(x0, x1, color=color, alpha=0.2, zorder=0)

        for x0 in zero_crossing_distances_km:
            for ax in axes:
                ax.axvline(
                    x0, color="gray", linestyle=":", linewidth=1, zorder=1
                )

        # Shade contiguous along-transect runs falling inside the
        # Coulomb-fit region (fit_mask), on every panel.
        fit_region_used = False
        fit_region_ranges = []
        fit_along_transect = None
        if fit_mask is not None:
            fit_along_transect = fit_mask[cell_indices]
            fit_change_indices = np.flatnonzero(
                np.diff(fit_along_transect.astype(np.int8)) != 0
            ) + 1
            fit_run_starts = np.concatenate(([0], fit_change_indices))
            fit_run_ends = np.concatenate(
                (fit_change_indices, [len(fit_along_transect) - 1])
            )
            for run_start, run_end in zip(fit_run_starts, fit_run_ends):
                if not fit_along_transect[run_start]:
                    continue
                fit_region_used = True
                x0 = distance_km[run_start]
                x1 = distance_km[run_end]
                fit_region_ranges.append((x0, x1))
                for ax in axes:
                    ax.axvspan(
                        x0, x1, color="tab:orange", alpha=0.2, zorder=0
                    )

        for ax, (label, field_spec) in zip(
            axes, all_fields.items()
        ):
            values, units, long_name = field_spec[:3]
            ylim = field_spec[3] if len(field_spec) > 3 else None
            if isinstance(values, dict):
                # Multi-line panel (e.g. bed & surface elevation).
                for line_label, arr in values.items():
                    ax.plot(
                        distance_km, arr[cell_indices],
                        label=line_label, zorder=2,
                    )
                ax.legend(loc="best", fontsize=8)
            elif label == implied_c_label and fit_along_transect is not None:
                # Only meaningful inside the Coulomb-fit region;
                # break the line elsewhere.
                sampled = values[cell_indices].astype(np.float64)
                sampled = np.where(fit_along_transect, sampled, np.nan)
                ax.plot(distance_km, sampled, zorder=2)
                if fitted_c is not None:
                    for i, (x0, x1) in enumerate(fit_region_ranges):
                        ax.hlines(
                            fitted_c, x0, x1,
                            color="k", linestyle="--", linewidth=1,
                            zorder=3,
                            label=(
                                f"fitted scalar C ({fitted_c:.3g})"
                                if i == 0 else None
                            ),
                        )
                    if fit_region_ranges:
                        ax.legend(loc="best", fontsize=8)
            else:
                ax.plot(distance_km, values[cell_indices], zorder=2)
            ax.set_ylabel(f"{label} [{units}]")
            ax.set_title(long_name, fontsize=10)
            ax.grid(True, alpha=0.3, zorder=1)
            if ylim is not None:
                # Clamp the autoscaled range to at most `ylim`, rather
                # than forcing it outright: this still prevents a
                # runaway/outlier-driven range (e.g. a fraction field
                # that should be in [0, 1]), while preserving natural
                # detail when the real data only span a narrow band
                # within those bounds (a hard set_ylim(*ylim) can
                # otherwise squash a genuinely narrow-but-valid range,
                # like 0.95-1.0, into an invisible sliver).
                ax.relim()
                ax.autoscale_view()
                data_lo, data_hi = ax.get_ylim()
                ax.set_ylim(
                    max(data_lo, ylim[0]), min(data_hi, ylim[1])
                )

            if (
                label == floatation_fraction_label
                and min_fraction_overburden is not None
            ):
                bound = min_fraction_overburden
                ax.axhline(
                    bound,
                    color="k",
                    linestyle="--",
                    linewidth=1,
                    zorder=3,
                    label=(
                        "min-fraction-overburden bound "
                        f"({bound:.3g})"
                    ),
                )
                ax.legend(loc="best", fontsize=8)

        axes[-1].set_xlabel("Along-transect distance [km]")

        if classes_used:
            legend_handles = [
                mpatches.Patch(
                    color=shading[cls][1], alpha=0.2, label=shading[cls][0]
                )
                for cls in sorted(classes_used)
            ]
        else:
            legend_handles = []
        if fit_region_used:
            legend_handles.append(
                mpatches.Patch(
                    color="tab:orange", alpha=0.2,
                    label="Coulomb C-fit region",
                )
            )
        if legend_handles:
            fig.legend(
                handles=legend_handles,
                loc="upper center",
                ncol=len(legend_handles),
                bbox_to_anchor=(0.5, 1.02),
                frameon=False,
            )

        if classes_used or fit_region_used:
            fig.suptitle(f"{name} transect", y=1.08)
        else:
            fig.suptitle(f"{name} transect")
        fig.tight_layout()

        out_path = os.path.join(plot_dir, f"transect_{name}.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

        print(f"Wrote transect plot: {out_path}")


def plot_maps(mesh_ds, fields, plot_dir, transects=None):
    """
    Plot Antarctic-wide maps of the given cell-centered fields on the
    native MALI mesh, using the `mosaic` package (available in the
    e3sm-unified conda environment) to render unstructured-mesh
    polygons directly in the mesh's native planar (polar
    stereographic) coordinates -- no reprojection needed since MALI's
    Antarctic meshes are already planar.

    Parameters
    ----------
    mesh_ds : xarray.Dataset
        A MALI mesh dataset containing the coordinate/connectivity
        arrays mosaic.Descriptor needs (xCell/yCell, verticesOnCell,
        cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnVertex) and
        the `on_a_sphere`/`is_periodic` global attributes. Antarctic
        MALI meshes are planar (on_a_sphere = "NO"), so no projection/
        transform is needed.
    fields : dict of str -> list of dict
        Mapping of output filename stem -> a list of one or more
        panels to plot side by side in that file. Each panel is a
        dict with keys:
            values : 1-D array, cell-centered values to plot.
            units : str, colorbar label.
            title : str, panel title.
            log : bool, optional (default False). Use a log-scale
                colormap (values <= 0 are masked out).
            cmap : str, optional (default "viridis"). Matplotlib
                colormap name.
            vmin, vmax : float, optional. Clip the colormap to this
                range (e.g. [0, 1] for a fraction field with a few
                out-of-range outliers, or [0.1, 10] for a log-scale
                ratio field diverging about 1.0). Also honored when
                `log` is True (as the LogNorm's vmin/vmax), not just
                for a linear-scale Normalize.
            mask : 1-D bool array, optional. Cells where mask is
                False are set to NaN (not displayed), e.g. to hide
                non-grounded cells on a grounded-ice-only field.
        e.g. {"bedRoughnessRC": [{"values": lam,
        "units": "Pa (m/yr)^-1/3",
        "title": "RC bed roughness Lambda", "log": True,
        "cmap": "turbo"}]}.
    plot_dir : str
        Directory to write output PNGs to (created if needed).
    transects : list of (str, ndarray, ndarray), optional
        (name, x, y) tuples, already projected into the mesh's planar
        CRS (e.g. via project_transect()), overlaid as lines on every
        panel of every map, with each transect's name labeled at its
        starting point. Omit/None to disable (default: no overlay).
    """
    try:
        import mosaic
    except ImportError as e:
        raise ImportError(
            "Plotting maps requires the 'mosaic' package. This is "
            "available in the e3sm-unified conda environment (see "
            "load_latest_e3sm_unified_*.sh); it is not required for "
            "the rest of this script."
        ) from e

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.colors
    import matplotlib.pyplot as plt

    os.makedirs(plot_dir, exist_ok=True)

    descriptor = mosaic.Descriptor(mesh_ds)

    for stem, panels in fields.items():
        fig, axes = plt.subplots(
            1, len(panels), figsize=(7 * len(panels), 7), squeeze=False
        )
        axes = axes[0]

        for ax, panel in zip(axes, panels):
            units = panel["units"]
            title = panel["title"]
            log = panel.get("log", False)
            cmap = panel.get("cmap", "viridis")
            vmin = panel.get("vmin")
            vmax = panel.get("vmax")
            mask = panel.get("mask")

            plot_values = np.asarray(panel["values"], dtype=np.float64)
            if mask is not None:
                plot_values = np.where(mask, plot_values, np.nan)

            norm = None
            if log:
                # Log-scale colormaps can't handle non-positive
                # values; mask them out (e.g. Lambda's full-Coulomb
                # reference cells, typically 0.0).
                plot_values = np.where(
                    plot_values > 0.0, plot_values, np.nan
                )
                norm = matplotlib.colors.LogNorm(
                    vmin=vmin if vmin is not None else np.nanmin(plot_values),
                    vmax=vmax if vmax is not None else np.nanmax(plot_values),
                )
            elif vmin is not None or vmax is not None:
                norm = matplotlib.colors.Normalize(
                    vmin=vmin, vmax=vmax, clip=True
                )
            array = xr.DataArray(plot_values, dims=("nCells",))
            coll = mosaic.polypcolor(
                ax, descriptor, array, cmap=cmap, norm=norm, aa=False
            )
            ax.set_aspect("equal")
            ax.set_title(title, fontsize=11)
            fig.colorbar(coll, ax=ax, label=units, shrink=0.8)

            if transects:
                for name, x, y in transects:
                    ax.plot(x, y, color="red", linewidth=1, zorder=3)
                    ax.annotate(
                        name, (x[0], y[0]), color="red", fontsize=7,
                        zorder=4,
                    )

        fig.tight_layout()

        out_path = os.path.join(plot_dir, f"map_{stem}.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

        print(f"Wrote map plot: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert MALI Weertman/Budd friction IC to Regularized Coulomb.",
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
        "--method",
        choices=["stress-match-fit", "transition-velocity"],
        default="stress-match-fit",
        help=(
            "How to derive bedRoughnessRC (Lambda) and the Regularized "
            "Coulomb friction coefficient (--mu-field) (default: "
            "stress-match-fit). \"stress-match-fit\" (the original "
            "method) fits a single, spatially-uniform scalar C over "
            "the fast-flowing region (speed > --critical-velocity), "
            "assumed already in the full-Coulomb regime, then solves "
            "per-cell for Lambda elsewhere so the RC law reproduces "
            "the source law's basal shear stress at each cell's "
            "current sliding speed (see fit_coulomb_C_fast_region()/ "
            "the Lambda solve below --critical-velocity/"
            "--lambda-reference-value are used by this method). "
            "\"transition-velocity\" instead picks a fixed transition "
            "velocity u0 (--transition-velocity) and computes, in "
            "closed form for every grounded cell, "
            "Lambda = u0 / (SECONDS_PER_YEAR * A * N^n) and a "
            "spatially-varying C = tau_b * (ub + u0)^(1/3) / "
            "(N * ub^(1/3)) (ub = current sliding speed, tau_b = "
            "source law shear stress mu*N_source*ub^qW), so that "
            "the RC law exactly reproduces the source law's stress "
            "at every grounded cell's current speed, with no "
            "fast/slow-region split (--transition-velocity/"
            "--mu-reference-value/"
            "--lambda-reference-value are used by this method)."
        )
    )
    parser.add_argument(
        "--critical-velocity", "--uc",
        type=float,
        default=None,
        help=(
            "Critical velocity u_c, e.g. in m/yr. Required, and only "
            "used, when --method=stress-match-fit."
        )
    )
    parser.add_argument(
        "--transition-velocity", "--u0",
        dest="transition_velocity",
        type=float,
        default=None,
        help=(
            "Transition velocity u0, e.g. in m/yr, used to derive "
            "Lambda = u0 / (SECONDS_PER_YEAR * A * N^n) and the "
            "per-cell Regularized Coulomb coefficient C = tau_b * "
            "(ub + u0)^(1/3) / (N * ub^(1/3)). Required, and only "
            "used, when --method=transition-velocity."
        )
    )
    parser.add_argument(
        "--mu-reference-value",
        type=float,
        default=0.3,
        help=(
            "Value assigned to the output --mu-field (the per-cell "
            "Regularized Coulomb coefficient C) at cells where the "
            "--method=transition-velocity closed-form solve is "
            "undefined (non-grounded, ice-free, or zero/invalid "
            "current sliding speed, N, or input mu -- the same cells "
            "that fall back to --lambda-reference-value for Lambda). "
            "Only used when --method=transition-velocity (default: "
            "0.3)."
        )
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

    parser.add_argument(
        "--lambda-reference-value",
        type=float,
        default=0.0,
        help=(
            "Value assigned to bedRoughnessRC (Lambda) at cells where "
            "no valid Lambda can be computed. For "
            "--method=stress-match-fit, this is all fast-flowing "
            "cells (speed > critical velocity, the same region used "
            "to fit C), plus any other grounded cell where the exact "
            "per-cell Lambda solve is ill-defined (source Tau_b "
            "already meets or exceeds the Coulomb limit C*N); Lambda "
            "-> 0 exactly reproduces the full-Coulomb limit there, so "
            "0.0 (default) is physically correct, though a small "
            "positive value can be used instead if a strictly-zero "
            "bed roughness is undesirable for other reasons. For "
            "--method=transition-velocity, this is any non-grounded, "
            "ice-free, or zero/invalid current-sliding-speed/N cell "
            "(default: 0.0)."
        )
    )

    extrapolate_terminus_group = parser.add_mutually_exclusive_group()
    extrapolate_terminus_group.add_argument(
        "--extrapolate-terminus-cells",
        dest="extrapolate_terminus_cells",
        action="store_true",
        default=False,
        help=(
            "Discard the computed bedRoughnessRC (Lambda) value (and, "
            "for --method=transition-velocity, the computed --mu-"
            "field/C value) at every grounded cell adjacent to the "
            "grounding line (a floating-ice neighbor) or to a "
            "grounded marine terminus (an ice-free, bed-below-sea-"
            "level neighbor, i.e. a tidewater-glacier-style calving "
            "front with no floating shelf), and at every non-grounded "
            "cell (floating ice, ice-free ocean, ice-free land) -- "
            "i.e. the entire mesh domain outside the grounded "
            "interior -- then refill all of those cells by creep-"
            "fill extrapolation (--creep-fill-method), sourced "
            "purely from the remaining grounded-interior cells -- "
            "see creep_fill_extrapolate(). The discarded marginal "
            "cells are often the least reliable (e.g. noisiest "
            "velocity/thickness/bed data, or most sensitive to the "
            "exact grounding-line position), so this discards them "
            "in favor of extrapolating from more interior, better-"
            "constrained cells, and additionally gives every non-"
            "grounded cell a physically-reasonable (rather than a "
            "flat reference) value. Cells filled this way are marked "
            "in the maskTerminusExtrapolated diagnostic field (see "
            "--diagnostics). Default: disabled (use the directly "
            "computed values everywhere)."
        )
    )
    extrapolate_terminus_group.add_argument(
        "--no-extrapolate-terminus-cells",
        dest="extrapolate_terminus_cells",
        action="store_false",
        help=(
            "Do not discard/extrapolate grounding-line/grounded-"
            "marine-terminus/non-grounded cells; use the directly "
            "computed values everywhere (default)."
        )
    )
    parser.add_argument(
        "--creep-fill-method",
        choices=["idw", "min"],
        default="idw",
        help=(
            "Extrapolation method used by --extrapolate-terminus-"
            "cells to fill discarded grounding-line/grounded-marine-"
            "terminus cells from neighboring valid cells (default: "
            "idw). \"idw\": inverse-distance-weighted average of "
            "valid neighbors. \"min\": minimum value among valid "
            "neighbors (matches MPAS-Tools' "
            "conversion_exodus_init_to_mpasli_mesh.py \"min\" "
            "extrapolation option). Only used when "
            "--extrapolate-terminus-cells is set."
        )
    )

    # Effective pressure N
    parser.add_argument(
        "--effective-pressure-type",
        choices=["downs-johnson", "ocean-connection", "transition"],
        default="downs-johnson",
        help=(
            "How to compute the effective pressure N (default: "
            "downs-johnson). \"downs-johnson\" reproduces Albany's "
            "own internal \"Hydrostatic Computed At Nodes\" Effective "
            "Pressure Type formula, with \"Use Pressurized Bed Above "
            "Sea Level: true\" (see --min-fraction-overburden/"
            "--pressure-length-scale). \"ocean-connection\" reproduces "
            "the same Albany \"Hydrostatic Computed At Nodes\" "
            "Effective Pressure Type but with \"Use Pressurized Bed "
            "Above Sea Level: false\" -- i.e. full ocean (hydrostatic) "
            "water pressure wherever bed is below sea level, with no "
            "inland tapering/floor, and no dependence on "
            "--min-fraction-overburden/--pressure-length-scale (see "
            "ocean_connection_effective_pressure()). Both "
            "\"downs-johnson\" and \"ocean-connection\" are computed "
            "by Albany itself at runtime from thickness/bed, so no N "
            "field needs to be supplied for either. \"transition\" "
            "computes N here using a near-ocean/inland-transition "
            "parameterization (see --min-fraction-overburden/"
            "--pressure-length-scale/--transition-h-ocean) and writes "
            "it to the output for reference; NOTE: Albany does not "
            "yet have a way to consume this precomputed N directly "
            "for the Regularized Coulomb law, so \"transition\" cannot "
            "currently be used to actually run Albany (offline "
            "evaluation only, pending upstream Albany support)."
        )
    )

    # Downs & Johnson / Albany parameters (required only when
    # --effective-pressure-type=downs-johnson).
    # Effective-pressure parameters shared between the "downs-johnson"
    # and "transition" parameterizations (see --effective-pressure-
    # type below): both use the same inland-floatation-fraction-floor
    # convention (min_fraction_overburden) and a length scale in
    # meters (pressure_length_scale), even though each
    # parameterization uses them in a different formula. Overloading
    # the same CLI flags across both types makes it convenient to
    # switch between them.
    parser.add_argument(
        "--min-fraction-overburden",
        type=float,
        default=None,
        help=(
            "Prescribed inland floatation-fraction floor. Matches "
            'Albany\'s "Minimum Fraction Overburden Pressure" '
            "convention exactly: far inland, N/Pice -> "
            "1 - min_fraction_overburden, i.e. the floatation "
            "fraction Pw/Pice approaches min_fraction_overburden "
            "itself (its minimum value; near the grounding line the "
            "floatation fraction rises toward 1 regardless of this "
            "parameter). Used by both --effective-pressure-type="
            "downs-johnson and --effective-pressure-type=transition "
            "(see downs_johnson_effective_pressure()/"
            "effective_pressure4()). Required in either case."
        )
    )
    parser.add_argument(
        "--pressure-length-scale",
        type=float,
        default=None,
        help=(
            "Length scale [m, same units as bedTopography] used by "
            "the effective-pressure parameterization selected via "
            "--effective-pressure-type. For downs-johnson, this is "
            'Albany\'s "Length Scale Factor" (width of the bed-'
            "elevation sigmoid). For transition, this is the "
            "distance over which half the remaining difference "
            "between the near-ocean value and the inland target "
            "fraction is removed moving inland (see "
            "effective_pressure4()); 0 disables the smooth "
            "transition (step function) in that case. The two "
            "parameterizations use this value in different formulas "
            "-- a value tuned for one is not automatically "
            "appropriate for the other. Required in either case."
        )
    )

    # "transition"-only effective-pressure parameter (required only
    # when --effective-pressure-type=transition).
    parser.add_argument(
        "--transition-h-ocean",
        type=float,
        default=25.0,
        help=(
            "Height above flotation [m, same units as bedTopography] "
            "below which N is assumed set purely by the "
            "ocean-connected fraction, with no inland transition "
            "applied (see effective_pressure4()); only used when "
            "--effective-pressure-type=transition (default: 25.0)"
        )
    )

    # -------------------------------------------------------------
    # Input source friction law: effective pressure N_source
    # -------------------------------------------------------------
    # These options are independent of --effective-pressure-type
    # above -- they describe the effective pressure implied by the
    # *input* MALI friction law (mu, muFriction), not Albany's own N
    # for the Regularized Coulomb law. See "Source friction law"
    # discussion in the module docstring.
    parser.add_argument(
        "--source-effective-pressure-type",
        choices=[
            "constant", "downs-johnson", "ocean-connection",
            "transition", "field",
        ],
        default="constant",
        help=(
            "How to compute the effective pressure N_source implied "
            "by the input source friction law "
            "(Tau_b_source = mu * N_source * speed^qW; default: "
            "constant). \"constant\" uses a single spatially-uniform "
            "value given by --source-effective-pressure (default "
            "1.0, which reproduces the classic MALI Weertman law "
            "exactly, since Tau_b_source then reduces to "
            "mu * speed^qW). \"downs-johnson\"/\"ocean-connection\"/"
            "\"transition\" compute N_source here using the same "
            "formulas available for --effective-pressure-type (see "
            "--source-min-fraction-overburden/"
            "--source-pressure-length-scale/"
            "--source-transition-h-ocean), letting the source law be "
            "a genuine Budd-type law (effective pressure specified "
            "rather than assumed 1) with its own, independently "
            "configured effective pressure -- it need not match "
            "--effective-pressure-type/its parameters. \"field\" "
            "reads a precomputed N_source field directly from the "
            "input file (see --source-effective-pressure-field), in "
            "physical units of Pa. In all non-constant cases, "
            "N_source is converted to the same kPa-equivalent scale "
            "as N_albany before use (see "
            "ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT)."
        )
    )
    parser.add_argument(
        "--source-effective-pressure",
        type=float,
        default=1.0,
        help=(
            "Spatially-uniform value of N_source to use when "
            "--source-effective-pressure-type=constant (default: "
            "1.0, the exact Weertman-equivalent value: "
            "Tau_b_source = mu * speed^qW). Used as-is with no unit "
            "conversion; only meaningful relative to --weertman-q "
            "and mu's own units."
        )
    )
    parser.add_argument(
        "--source-effective-pressure-field",
        default=None,
        help=(
            "Name of a variable in the input file holding a "
            "precomputed N_source field [Pa], used only when "
            "--source-effective-pressure-type=field."
        )
    )
    parser.add_argument(
        "--source-min-fraction-overburden",
        type=float,
        default=None,
        help=(
            "Same convention as --min-fraction-overburden, but for "
            "the source law's own N_source. Used by both "
            "--source-effective-pressure-type=downs-johnson and "
            "--source-effective-pressure-type=transition. Required "
            "in either case; independent of --min-fraction-"
            "overburden (may use a different value)."
        )
    )
    parser.add_argument(
        "--source-pressure-length-scale",
        type=float,
        default=None,
        help=(
            "Same convention as --pressure-length-scale, but for the "
            "source law's own N_source. Used by both "
            "--source-effective-pressure-type=downs-johnson and "
            "--source-effective-pressure-type=transition. Required "
            "in either case; independent of --pressure-length-scale "
            "(may use a different value)."
        )
    )
    parser.add_argument(
        "--source-transition-h-ocean",
        type=float,
        default=25.0,
        help=(
            "Same convention as --transition-h-ocean, but for the "
            "source law's own N_source; only used when "
            "--source-effective-pressure-type=transition (default: "
            "25.0)."
        )
    )

    parser.add_argument("--rho-ice", type=float, default=910.0)
    parser.add_argument(
        "--rho-water",
        type=float,
        default=1028.0,
        help=(
            "Ocean/seawater density [kg m^-3] (default: 1028.0), "
            "used for ocean-connectivity terms: the grounded-ice "
            "test, the marine term in the effective-pressure "
            "parameterizations, and the floatation-fraction "
            "denominator. Distinct from --rho-freshwater, used for "
            "the subglacial water column in the hydropotential field."
        )
    )
    parser.add_argument(
        "--rho-freshwater",
        type=float,
        default=1000.0,
        help=(
            "Freshwater density [kg m^-3] (default: 1000.0), used "
            "only for the bed-elevation term of the Shreve hydraulic "
            "potential (hydropotential = rho_freshwater * gravity * "
            "bedTopography + Pw), matching the standard convention "
            "that subglacial water is fresh, distinct from the "
            "ocean/seawater density (--rho-water) used elsewhere."
        )
    )
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
        help=(
            "Source-law friction field (default: muFriction). Its "
            "correct physical units -- per a correction to "
            "MPAS-Tools' Registry.xml -- are kPa * (yr/m)^qW when "
            "--source-effective-pressure-type=constant (not "
            "Pa * (yr/m)^qW), or plain (yr/m)^qW otherwise (see "
            "--source-effective-pressure-type); this script's use of "
            "mu is dimensionally consistent with that (see "
            "fit_coulomb_C_fast_region())."
        )
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
        "--ascii-mesh-dir",
        default=".",
        help=(
            "Directory containing MALI's `config_write_albany_ascii_"
            "mesh` ascii output (`thickness.ascii`, "
            "`bed_topography.ascii`, `mpas_cellID.ascii`; default: "
            "'.'). MALI extends the input mesh by one cell around "
            "its boundary when building the FEM mesh handed to "
            "Albany, and may assign a fixed minimum thickness to "
            "that extended ring, so these files -- not "
            "--thickness-field/--bed-field -- are always used as the "
            "authoritative ice thickness/bed elevation for every "
            "downstream calculation (effective pressure, grounded/ "
            "floating classification, grounding-line/terminus "
            "identification, diagnostics, and plots)."
        )
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
        "--floatation-fraction-field",
        default="floatationFraction",
        help="Name for diagnostic floatation fraction field (Pw/Pice)"
    )
    parser.add_argument(
        "--hydropotential-field",
        default="hydropotential",
        help="Name for diagnostic hydraulic potential field"
    )

    parser.add_argument(
        "--transects-dir",
        default=None,
        help=(
            "Path to a geometric_features landice/transect directory "
            "(e.g. .../geometric_features/geometric_data/landice/"
            "transect), containing one subdirectory per named "
            "transect, each with a transect.geojson (LineString "
            "lon/lat, degrees). Required if --plot-transects is set."
        )
    )
    parser.add_argument(
        "--transect-names",
        nargs="+",
        default=[
            "Thwaites", "Totten", "Jutulstraumen", "Foundation",
            "Bindschadler", "Pine_Island",
        ],
        help=(
            "Names of transects to plot (subdirectory names under "
            "--transects-dir). Default: Thwaites Totten Jutulstraumen "
            "Foundation Bindschadler Pine_Island."
        )
    )
    parser.add_argument(
        "--plot-dir",
        default="diagnostic_plots",
        help=(
            "Directory to write all diagnostic PNG plots to (transect "
            "plots and Antarctic-wide maps, all directly in this "
            "single directory, no subdirectories). Default: "
            "diagnostic_plots, created if needed."
        )
    )

    plot_transects_group = parser.add_mutually_exclusive_group()
    plot_transects_group.add_argument(
        "--plot-transects",
        dest="plot_transects",
        action="store_true",
        default=False,
        help=(
            "Plot effectivePressure, floatationFraction, and "
            "hydropotential along the transects in --transect-names, "
            "sampled at MALI cell centers nearest each transect "
            "point (projected from lon/lat into the MALI mesh's "
            "planar CRS, EPSG:3031 for Antarctica). Requires "
            "--diagnostics (the fields plotted are diagnostic "
            "fields) and --transects-dir, plus the pyproj/scipy/"
            "matplotlib packages. Also automatically produces "
            "Antarctic-wide maps (in --plot-dir) of the same "
            "diagnostic fields, the diagnostic masks, bedRoughnessRC "
            "(log scale), and the original muFriction field (log "
            "scale, for comparison), using the 'mosaic' package "
            "(available in the e3sm-unified conda environment). "
            "Default: disabled."
        )
    )
    plot_transects_group.add_argument(
        "--no-plot-transects",
        dest="plot_transects",
        action="store_false",
        help="Do not plot transects (default)."
    )

    plot_transects_on_maps_group = parser.add_mutually_exclusive_group()
    plot_transects_on_maps_group.add_argument(
        "--plot-transects-on-maps",
        dest="plot_transects_on_maps",
        action="store_true",
        default=True,
        help=(
            "Overlay the --transect-names transect lines (with "
            "labels) on every panel of every Antarctic-wide map "
            "produced by --plot-transects. Default: enabled."
        )
    )
    plot_transects_on_maps_group.add_argument(
        "--no-plot-transects-on-maps",
        dest="plot_transects_on_maps",
        action="store_false",
        help="Do not overlay transect lines on the Antarctic-wide maps."
    )

    diagnostics_group = parser.add_mutually_exclusive_group()
    diagnostics_group.add_argument(
        "--diagnostics",
        dest="diagnostics",
        action="store_true",
        default=True,
        help=(
            "Include diagnostic fields (effectivePressure, flowRateA, "
            "floatationFraction, hydropotential, maskGrounded, "
            "maskFastFlowing, maskValidBedRoughnessRC) in the output "
            "NetCDF, useful for debugging/evaluation/visualization "
            "(default: enabled)."
        )
    )
    diagnostics_group.add_argument(
        "--no-diagnostics",
        dest="diagnostics",
        action="store_false",
        help=(
            "Omit diagnostic fields from the output NetCDF, e.g. for "
            "production runs where only bedRoughnessRC and muFriction "
            "are needed and extra fields are unwanted clutter."
        )
    )

    parser.add_argument(
        "--time-index",
        type=int,
        default=0,
        help="Time index for Time-dependent IC fields (default: 0)"
    )

    args = parser.parse_args()

    if args.method == "stress-match-fit":
        if args.critical_velocity is None:
            parser.error(
                "--critical-velocity is required when "
                "--method=stress-match-fit"
            )
        if args.transition_velocity is not None:
            parser.error(
                "--transition-velocity is only used with "
                "--method=transition-velocity (got "
                "--method=stress-match-fit)"
            )
    else:
        if args.transition_velocity is None:
            parser.error(
                "--transition-velocity is required when "
                "--method=transition-velocity"
            )
        if args.critical_velocity is not None:
            parser.error(
                "--critical-velocity is only used with "
                "--method=stress-match-fit (got "
                "--method=transition-velocity)"
            )

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

    if args.effective_pressure_type != "ocean-connection":
        if (
            args.min_fraction_overburden is None
            or args.pressure_length_scale is None
        ):
            parser.error(
                "--min-fraction-overburden and --pressure-length-scale "
                "are required (used by both "
                "--effective-pressure-type=downs-johnson and "
                "--effective-pressure-type=transition; not needed for "
                "--effective-pressure-type=ocean-connection)"
            )

    if args.source_effective_pressure_type == "field":
        if args.source_effective_pressure_field is None:
            parser.error(
                "--source-effective-pressure-field is required when "
                "--source-effective-pressure-type=field"
            )
    elif args.source_effective_pressure_field is not None:
        parser.error(
            "--source-effective-pressure-field is only used with "
            "--source-effective-pressure-type=field"
        )

    if args.source_effective_pressure_type in ("downs-johnson", "transition"):
        if (
            args.source_min_fraction_overburden is None
            or args.source_pressure_length_scale is None
        ):
            parser.error(
                "--source-min-fraction-overburden and "
                "--source-pressure-length-scale are required (used by "
                "both --source-effective-pressure-type=downs-johnson "
                "and --source-effective-pressure-type=transition)"
            )
    elif (
        args.source_min_fraction_overburden is not None
        or args.source_pressure_length_scale is not None
    ):
        parser.error(
            "--source-min-fraction-overburden and "
            "--source-pressure-length-scale are only used with "
            "--source-effective-pressure-type=downs-johnson or "
            "=transition"
        )

    if args.plot_transects:
        if not args.diagnostics:
            parser.error(
                "--plot-transects requires --diagnostics (the fields "
                "plotted -- effectivePressure, floatationFraction, "
                "hydropotential -- are only computed when diagnostics "
                "are enabled)"
            )
        if args.transects_dir is None:
            parser.error(
                "--transects-dir is required when --plot-transects "
                "is set"
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
    if args.extrapolate_terminus_cells:
        # MPAS mesh connectivity, needed to identify grounding-line/
        # grounded-marine-terminus cells and to creep-fill them.
        required.extend(["cellsOnCell", "nEdgesOnCell", "xCell", "yCell"])

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
    # Override thickness/bed with the values Albany actually used,
    # from MALI's `config_write_albany_ascii_mesh` ascii output
    # (--ascii-mesh-dir). MALI extends the input mesh by one cell
    # around its boundary when building the FEM mesh for Albany
    # (assigning a fixed minimum thickness to that extended ring), so
    # these ascii-derived values -- not the plain netCDF
    # thickness/bedTopography -- are used for every downstream
    # calculation. Cells outside the ascii files' coverage (i.e. not
    # part of the FEM/Albany domain) keep their netCDF values.
    # -------------------------------------------------------------
    ascii_cell_index, ascii_H, ascii_bed = load_albany_ascii_geometry(
        args.ascii_mesh_dir, n_cells=H.shape[0]
    )
    n_overridden = ascii_cell_index.shape[0]
    n_changed = int(np.count_nonzero(
        (H[ascii_cell_index] != ascii_H)
        | (bed[ascii_cell_index] != ascii_bed)
    ))
    H[ascii_cell_index] = ascii_H
    bed[ascii_cell_index] = ascii_bed
    print(
        f"Overrode thickness/bedTopography for {n_overridden} cells "
        f"from Albany ascii mesh files in {args.ascii_mesh_dir!r} "
        f"({n_changed} differ from the input netCDF file, e.g. due "
        "to MALI's one-cell mesh extension/fixed margin thickness)"
    )

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
    if args.effective_pressure_type == "downs-johnson":
        N = downs_johnson_effective_pressure(
            thickness=H,
            bed=bed,
            min_fraction_overburden=args.min_fraction_overburden,
            length_scale=args.pressure_length_scale,
            rho_i=args.rho_ice,
            rho_w=args.rho_water,
            gravity=args.gravity,
        )
    elif args.effective_pressure_type == "ocean-connection":
        N = ocean_connection_effective_pressure(
            thickness=H,
            bed=bed,
            rho_i=args.rho_ice,
            rho_w=args.rho_water,
            gravity=args.gravity,
        )
    else:
        N = effective_pressure4(
            thickness=H,
            bed=bed,
            min_fraction_overburden=args.min_fraction_overburden,
            length_scale=args.pressure_length_scale,
            rho_i=args.rho_ice,
            rho_w=args.rho_water,
            gravity=args.gravity,
            h_ocean=args.transition_h_ocean,
        )

    # Albany's own internal effective pressure (e.g. "Hydrostatic"
    # Effective Pressure Type) is computed from bed/thickness fields
    # that MALI's coupling interface has already divided by 1000 (m ->
    # km) before Albany sees them, which numerically makes Albany's
    # internal N equal to physical N[Pa] / 1000. C must be derived
    # against this same scale for consistency with Albany's own
    # "beta = C * N * |u|^(q-1)" evaluation.
    N_albany = N / ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT

    # -------------------------------------------------------------
    # Effective pressure implied by the *input source* friction law,
    # N_source (see --source-effective-pressure-type). Independent of
    # N/N_albany above -- may use a different parameterization and/or
    # different parameter values. Always expressed in the same
    # kPa-equivalent scale as N_albany (physical Pa / 1000) so that
    # Tau_b_source = mu * N_source * speed^qW comes out on the same
    # scale as N_albany with no further conversion, except for
    # "constant", whose value is used as-is (default 1.0 exactly
    # reproduces classic Weertman).
    # -------------------------------------------------------------
    if args.source_effective_pressure_type == "constant":
        N_source_kpa = np.full_like(H, args.source_effective_pressure)
    elif args.source_effective_pressure_type == "field":
        N_source_kpa = (
            cell_field(args.source_effective_pressure_field)
            / ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT
        )
    else:
        N_source_kpa = compute_named_effective_pressure(
            args.source_effective_pressure_type,
            thickness=H,
            bed=bed,
            rho_i=args.rho_ice,
            rho_w=args.rho_water,
            gravity=args.gravity,
            min_fraction_overburden=args.source_min_fraction_overburden,
            pressure_length_scale=args.source_pressure_length_scale,
            transition_h_ocean=args.source_transition_h_ocean,
        ) / ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT

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
    # Grounding-line / grounded-marine-terminus cell identification
    # (--extrapolate-terminus-cells only): grounded cells immediately
    # adjacent to floating ice (the grounding line proper) or to
    # ice-free, bed-below-sea-level ocean (a grounded ice front with
    # no floating shelf, e.g. a tidewater glacier calving front).
    # Land-terminating margins (ice-free, bed at/above sea level) are
    # deliberately not included -- these are neither a grounding line
    # nor a marine terminus.
    #
    # These cells, plus every non-grounded cell (floating ice,
    # ice-free ocean, ice-free land), make up `fill_mask`: the
    # portion of the *entire mesh domain* whose values are discarded
    # and creep-filled by extrapolation, sourced from `keep_mask`
    # (the grounded interior, i.e. grounded ice minus its outermost
    # terminus row).
    # -------------------------------------------------------------
    if args.extrapolate_terminus_cells:
        ice_free = H <= 0.0
        floating = (~ice_free) & (~grounded)
        ice_free_ocean = ice_free & (bed < 0.0)

        cells_on_cell = np.asarray(ds["cellsOnCell"].values)
        n_edges_on_cell = np.asarray(ds["nEdgesOnCell"].values)
        x_cell = np.asarray(ds["xCell"].values, dtype=np.float64)
        y_cell = np.asarray(ds["yCell"].values, dtype=np.float64)

        terminus_cells = grounded & cell_has_neighbor_where(
            floating | ice_free_ocean, cells_on_cell, n_edges_on_cell
        )
        keep_mask = grounded & ~terminus_cells
        fill_mask = ~keep_mask

        print(
            "Grounding-line/grounded-marine-terminus cells discarded "
            f": {int(np.count_nonzero(terminus_cells))} "
            f"/ {int(np.count_nonzero(grounded))} grounded"
        )
        print(
            "Total cells to be extrapolated over (terminus + all "
            f"non-grounded cells) : {int(np.count_nonzero(fill_mask))} "
            f"/ {H.shape[0]} total cells"
        )
    else:
        terminus_cells = None
        keep_mask = None
        fill_mask = None

    # -------------------------------------------------------------
    # Basal shear stress implied by the input source friction law at
    # each cell's current sliding speed (from --velocity-x-field/
    # --velocity-y-field). Used by both --method options below to
    # solve for Lambda/C so that the RC law reproduces this same
    # stress:
    #
    #   Source law:          beta_source = mu * N_source * u^(qW-1)
    #                        Tau_b_source = beta_source * u
    #                                     = mu * N_source * u^qW
    #
    # (N_source == 1 everywhere is the classic Weertman law; see
    # --source-effective-pressure-type.)
    # -------------------------------------------------------------
    speed_defined = (
        grounded
        & np.isfinite(N) & (N > 0.0)
        & np.isfinite(mu) & (mu > 0.0)
        & np.isfinite(N_source_kpa) & (N_source_kpa > 0.0)
        & np.isfinite(speed) & (speed > 0.0)
    )

    tau_b_source = np.full_like(N, np.nan)
    tau_b_source[speed_defined] = (
        mu[speed_defined] * N_source_kpa[speed_defined]
        * speed[speed_defined] ** args.weertman_q
    )

    if args.method == "stress-match-fit":
        # ---------------------------------------------------------
        # Optimal C: fit against the fast-flowing region only,
        # assuming it is already in the fully-plastic Coulomb regime
        # of the RC law (Tau_b = C * N).
        # ---------------------------------------------------------
        fast_flowing = grounded & (speed > args.critical_velocity)

        C, fit_mask = fit_coulomb_C_fast_region(
            tau_b_source=tau_b_source,
            N=N_albany,
            area=area,
            mask=fast_flowing,
        )

        # ---------------------------------------------------------
        # Lambda / Albany Bed Roughness
        #
        # Rather than assuming the sliding speed equals a prescribed
        # critical velocity, Lambda is solved for exactly, per cell,
        # by requiring that Albany's Regularized Coulomb law
        # reproduce the same Tau_b computed above.
        #
        #   Regularized Coulomb (LandIce_BasalFrictionCoefficient_Def.hpp):
        #                        beta_RC = C * N * u^(p-1)
        #                                  / (u + Lambda*scaling*A*N^n)^p
        #                        Tau_b   = beta_RC * u
        #                                = C * N * u^p
        #                                  / (u + Lambda*scaling*A*N^n)^p
        #
        #   Setting Tau_b_RC == Tau_b_source and solving for Lambda:
        #
        #     u + Lambda*scaling*A*N^n = u * (C*N / Tau_b_source)^(1/p)
        #                    = u * (C*N / (mu * N_source * u^qW))^(1/p)
        #
        #     Lambda[m] = u * [(C*N / (mu * N_source * u^qW))^(1/p) - 1]
        #                 / (SECONDS_PER_YEAR * A * N^n)
        #     Lambda[m] is the value actually stored in the output field
        #     -- see the module docstring / comment above
        #     SECONDS_PER_YEAR's definition for why no further m -> km
        #     conversion is applied here (MALI's own coupling interface
        #     performs that conversion before Albany ever sees the field).
        #
        # where u is in m/yr, A is in Pa^-3 s^-1, N (raw Pa) is used
        # here exactly as in the previous uc-based derivation
        # (Albany's internal km/kPa/yr "scaling" factor reduces to
        # the plain SECONDS_PER_YEAR factor once Lambda is expressed
        # in meters, N in Pa, and the meters-based Lambda is stored
        # as-is, in meters, with MALI's coupling interface performing
        # the m -> km conversion before Albany sees it). N_albany
        # (the kPa-equivalent convention) is used for the "C*N"
        # Coulomb-limit term, matching how C was itself fit.
        #
        # Because Albany's Regularized Coulomb law can never produce
        # a shear stress above the Coulomb limit C*N (attained only
        # in the Lambda -> 0 limit), cells where the source law's
        # Tau_b at the current speed already meets or exceeds C*N
        # have no valid (non-negative) solution for Lambda.
        # Physically, Lambda -> 0 is *exactly* the fully-plastic
        # Coulomb regime (Tau_b_RC saturates at its maximum
        # achievable value, C*N, independent of speed), so setting
        # Lambda = args.lambda_reference_value at these cells is not
        # an arbitrary filler value -- it is the correct behavior for
        # cells that the fast-flowing/full-Coulomb assumption is
        # designed to describe in the first place.
        #
        # Fast-flowing cells (the same region used to fit C, i.e.
        # speed > critical_velocity) are *always* assumed to be in
        # this full-Coulomb regime and are therefore always forced to
        # Lambda = args.lambda_reference_value, regardless of what
        # the per-cell algebraic solve above would otherwise give --
        # the exact per-cell solve is not attempted there at all,
        # since matching the source law exactly at high speed is
        # not the goal (the fast region is assumed C-limited by
        # construction).
        # ---------------------------------------------------------
        Lambda = np.full_like(N, args.lambda_reference_value)

        valid_speed = speed_defined & ~fast_flowing

        stress_ratio = np.full_like(N, np.nan)
        stress_ratio[valid_speed] = (
            (C * N_albany[valid_speed]) / tau_b_source[valid_speed]
        )

        lambda_mask = (
            valid_speed & np.isfinite(stress_ratio) & (stress_ratio > 1.0)
        )

        Lambda[lambda_mask] = (
            speed[lambda_mask]
            * (stress_ratio[lambda_mask] ** (1.0 / RC_POWER_EXPONENT) - 1.0)
            / (SECONDS_PER_YEAR * A[lambda_mask] * N[lambda_mask] ** args.glen_n)
        )

        n_unreachable = int(np.count_nonzero(valid_speed & ~lambda_mask))
        if n_unreachable > 0:
            print(
                f"NOTE: {n_unreachable} slow-flowing grounded cells have a "
                "source-law basal shear stress (mu*N_source*u^qW) at the "
                "current sliding speed that meets or exceeds the Coulomb "
                f"limit C*N; Lambda set to {args.lambda_reference_value:g} "
                "(maximal Coulomb sliding) at these cells."
            )
        print(
            f"Fast-flowing cells forced to Lambda = "
            f"{args.lambda_reference_value:g} (full-Coulomb assumption) : "
            f"{int(np.count_nonzero(fast_flowing))}"
        )

        # Floating/ice-free/invalid-speed/unreachable-stress cells
        # are deliberately left at zero (see Lambda initialization
        # above).

        # -----------------------------------------------------------
        # Diagnostics
        # -----------------------------------------------------------
        # "Local" implied C: the per-cell ratio of the actual
        # source-law shear stress (at the cell's current sliding
        # speed) to N, i.e. what C would have to be for that cell
        # alone to be exactly in the full-Coulomb regime. Its spread
        # within the fast-flowing fit region gives a sense of how
        # well a single scalar C fits that region. Only meaningful
        # (and only plotted) for --method=stress-match-fit, since
        # --method=transition-velocity already writes a spatially-
        # varying, exactly-matching C to --mu-field directly.
        local_C = np.full_like(N, np.nan)
        local_C[speed_defined] = (
            tau_b_source[speed_defined] / N_albany[speed_defined]
        )
    else:
        # ---------------------------------------------------------
        # Transition-velocity method: Lambda and a spatially-varying
        # C are both computed in closed form from a fixed transition
        # velocity u0, with no fast/slow-region split -- see
        # solve_transition_velocity().
        # ---------------------------------------------------------
        fast_flowing = grounded & (speed > args.transition_velocity)
        fit_mask = None
        local_C = None

        Lambda, C = solve_transition_velocity(
            tau_b_source=tau_b_source,
            speed=speed,
            N=N,
            N_albany=N_albany,
            A=A,
            glen_n=args.glen_n,
            transition_velocity=args.transition_velocity,
            lambda_reference_value=args.lambda_reference_value,
            mu_reference_value=args.mu_reference_value,
            valid=speed_defined,
        )

        # Every grounded cell with a well-defined speed/N/mu gets an
        # exact closed-form solve (no fast/slow-region split, unlike
        # --method=stress-match-fit).
        lambda_mask = speed_defined

        print(
            f"Grounded cells with valid closed-form Lambda/C solve : "
            f"{int(np.count_nonzero(speed_defined))} "
            f"/ {int(np.count_nonzero(grounded))} grounded"
        )
        print(
            f"Cells above transition velocity, u0 (maskFastFlowing) : "
            f"{int(np.count_nonzero(fast_flowing))}"
        )

    # -------------------------------------------------------------
    # Whole-domain extrapolation (--extrapolate-terminus-cells):
    # discard the just-computed Lambda (and, for
    # --method=transition-velocity, C) everywhere outside the
    # grounded interior (`fill_mask` = grounding-line/grounded-
    # marine-terminus cells plus every non-grounded cell -- floating
    # ice, ice-free ocean, ice-free land) and creep-fill them from
    # `keep_mask` (the grounded interior) -- see
    # creep_fill_extrapolate().
    # -------------------------------------------------------------
    if args.extrapolate_terminus_cells and np.any(fill_mask):
        lambda_before = Lambda[fill_mask].copy()
        Lambda = creep_fill_extrapolate(
            Lambda,
            keep_mask=keep_mask,
            fill_mask=fill_mask,
            x_cell=x_cell,
            y_cell=y_cell,
            cells_on_cell=cells_on_cell,
            n_edges_on_cell=n_edges_on_cell,
            method=args.creep_fill_method,
        )
        print(
            "Whole-domain bedRoughnessRC (Lambda) discarded and "
            f"extrapolated ({args.creep_fill_method}) from the "
            "grounded interior: "
            f"{np.nanmin(lambda_before):.6e} -- "
            f"{np.nanmax(lambda_before):.6e} (before) -> "
            f"{np.nanmin(Lambda[fill_mask]):.6e} -- "
            f"{np.nanmax(Lambda[fill_mask]):.6e} (after)"
        )

        if args.method == "transition-velocity":
            c_before = C[fill_mask].copy()
            C = creep_fill_extrapolate(
                C,
                keep_mask=keep_mask,
                fill_mask=fill_mask,
                x_cell=x_cell,
                y_cell=y_cell,
                cells_on_cell=cells_on_cell,
                n_edges_on_cell=n_edges_on_cell,
                method=args.creep_fill_method,
            )
            print(
                "Whole-domain muFriction (C) discarded and "
                f"extrapolated ({args.creep_fill_method}) from the "
                "grounded interior: "
                f"{np.nanmin(c_before):.6e} -- "
                f"{np.nanmax(c_before):.6e} (before) -> "
                f"{np.nanmin(C[fill_mask]):.6e} -- "
                f"{np.nanmax(C[fill_mask]):.6e} (after)"
            )
        # Not applied to --method=stress-match-fit's C: that C is a
        # single scalar broadcast to every cell, so extrapolation
        # would have no effect.

        # Cells re-derived by extrapolation are no longer a "valid
        # solve" in the maskValidBedRoughnessRC sense (they were
        # deliberately discarded, not solved), but they are also not
        # simply left at a reference value -- track them separately
        # via maskTerminusExtrapolated (see the diagnostics section
        # below) rather than folding them into lambda_mask.

    # -------------------------------------------------------------
    # Regularized Coulomb regime ratio: u / (Lambda * A * N^n),
    # dimensionally consistent via the same SECONDS_PER_YEAR factor
    # used everywhere else in this script to relate Lambda (solved/
    # stored in meters) to a velocity scale (see the Lambda solves
    # above and solve_transition_velocity()). The denominator,
    # SECONDS_PER_YEAR * Lambda * A * N^n, is exactly the implied
    # critical/transition velocity u_c at which the Regularized
    # Coulomb law's basal shear stress switches over between its two
    # asymptotic regimes:
    #
    #   beta_RC = C * N * u^(p-1) / (u + u_c)^p
    #
    #   u >> u_c (ratio >> 1): (u + u_c)^p ~ u^p, so
    #       Tau_b = beta_RC * u -> C * N, independent of u -- the
    #       fully-plastic Coulomb regime.
    #   u << u_c (ratio << 1): (u + u_c)^p ~ u_c^p, so
    #       Tau_b ~ (C * N / u_c^p) * u^p -- a power-law (Weertman-
    #       like) regime.
    #
    # So this ratio is a direct, per-cell diagnostic of how close a
    # cell's current sliding speed places it to either regime:
    # ratio >> 1 is Coulomb-like, ratio << 1 is power-law-like,
    # ratio ~ 1 is the transition itself.
    # -------------------------------------------------------------
    with np.errstate(divide="ignore", invalid="ignore"):
        critical_velocity_implied = (
            SECONDS_PER_YEAR * Lambda * A * N ** args.glen_n
        )
        regime_ratio = speed / critical_velocity_implied

    print()
    print("MALI Weertman/Budd -> Regularized Coulomb conversion")
    print("------------------------------------------------")
    print(f"Input file                    : {args.input}")
    print(f"Method                        : {args.method}")
    if args.method == "stress-match-fit":
        print(f"Critical velocity, uc         : {args.critical_velocity:g}")
    else:
        print(f"Transition velocity, u0       : {args.transition_velocity:g}")
    print(f"Source power exponent, qW     : {args.weertman_q:g}")
    print(f"Source effective pressure type: {args.source_effective_pressure_type}")
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
    if args.method == "stress-match-fit":
        print(
            f"Fast-flowing (speed > uc) cells used in C fit : "
            f"{np.count_nonzero(fit_mask)} / {int(np.count_nonzero(grounded))} grounded"
        )
        print(f"Fast-flowing area used         : {np.sum(area[fit_mask]):.10e}")
        print()
        print(f"Optimal C                     : {C:.16e}")
        print()
        print(
            "N range on fit domain        : "
            f"{np.nanmin(N[fit_mask]):.6e} -- "
            f"{np.nanmax(N[fit_mask]):.6e}"
        )
    else:
        print()
        print(
            "C range (spatially varying)  : "
            + (
                f"{np.nanmin(C[speed_defined]):.6e} -- "
                f"{np.nanmax(C[speed_defined]):.6e}"
                if np.any(speed_defined) else "n/a (no valid cells)"
            )
        )
        print()
    print(
        "Basal sliding speed range (all grounded) : "
        + (
            f"{np.nanmin(speed[speed_defined]):.6e} -- "
            f"{np.nanmax(speed[speed_defined]):.6e} m/yr"
            if np.any(speed_defined) else "n/a (no valid cells)"
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
    if args.method == "stress-match-fit":
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

    # The source-law muFriction field is not used directly by the
    # Regularized Coulomb law; overwrite it with the RC coefficient C
    # (Albany's "Mu"/"Mu Field Name"), so that an Albany YAML using
    # "Mu Type: Field" (reading this same field name, per --mu-field)
    # picks up the fitted/calculated value. For
    # --method=stress-match-fit, C is a single scalar, area-weighted
    # fit broadcast uniformly to every cell (see
    # fit_coulomb_C_fast_region()); for --method=transition-velocity,
    # C already varies per cell (see solve_transition_velocity()).
    if args.method == "stress-match-fit":
        mu_field_values = np.full_like(N, C)
        mu_field_long_name = (
            "Albany regularized-Coulomb coefficient C (Mu), "
            "spatially uniform, area-weighted fit over the "
            "fast-flowing region -- see fit_coulomb_C_fast_region()"
        )
    else:
        mu_field_values = C
        mu_field_long_name = (
            "Albany regularized-Coulomb coefficient C (Mu), computed "
            "in closed form per cell from the transition velocity u0 "
            "-- see solve_transition_velocity()"
        )

    if args.extrapolate_terminus_cells and args.method == "transition-velocity":
        mu_field_long_name += (
            "; the grounding-line/grounded-marine-terminus band and "
            "every non-grounded cell (maskTerminusExtrapolated) were "
            "discarded and creep-fill extrapolated "
            f"({args.creep_fill_method}) from the grounded interior "
            "-- see --extrapolate-terminus-cells/"
            "creep_fill_extrapolate()"
        )

    out[args.mu_field] = xr.DataArray(
        mu_field_values,
        dims=(ncell_dim,),
        attrs={
            "long_name": mu_field_long_name,
            "units": "1",
        },
    )

    if args.method == "stress-match-fit":
        lambda_field_description = (
            "Fast-flowing cells (speed > critical velocity) and "
            "any other grounded cell where the solve below is "
            "ill-defined are assumed to be in the full-Coulomb "
            f"regime and set to {args.lambda_reference_value:g} "
            "(see maskFastFlowing/maskValidBedRoughnessRC). "
            "Elsewhere, Lambda is solved exactly so that the "
            "Regularized Coulomb law reproduces the source law's "
            "basal shear stress (mu*N_source*u^qW) at the cell's "
            "actual current sliding "
            "speed u (from velocity-x/y-field, last "
            "nVertInterfaces level): Lambda[m] = u * "
            "[(C*N/(mu*N_source*u^qW))^(1/qR) - 1] "
            "/ (SECONDS_PER_YEAR * A "
            "* N^n), with u in m/yr, A in Pa^-3 s^-1, N in Pa, "
            "matching Albany's internal secsInYr scaling in "
            "LandIce_BasalFrictionCoefficient_Def.hpp; the stored "
            "value is Lambda in meters, as-is (matching MALI's "
            "Registry.xml units=\"m\" declaration for this field) "
            "-- MALI's own Albany coupling interface "
            "(Interface_velocity_solver.cpp) divides this field "
            "by 1000 before Albany sees it, satisfying Albany's "
            "'bedRoughness in km' scaling convention"
        )
    else:
        lambda_field_description = (
            "Non-grounded, ice-free, or zero/invalid current-"
            f"sliding-speed/N/mu cells are set to "
            f"{args.lambda_reference_value:g} (see "
            "maskValidBedRoughnessRC). Elsewhere (every grounded "
            "cell with a valid solve), Lambda = u0 / "
            "(SECONDS_PER_YEAR * A * N^n), with u0 the transition "
            "velocity (--transition-velocity), A in Pa^-3 s^-1, N in "
            "Pa -- see solve_transition_velocity(); the stored value "
            "is Lambda in meters, as-is (matching MALI's Registry.xml "
            "units=\"m\" declaration for this field) -- MALI's own "
            "Albany coupling interface (Interface_velocity_solver.cpp) "
            "divides this field by 1000 before Albany sees it, "
            "satisfying Albany's 'bedRoughness in km' scaling "
            "convention"
        )

    if args.extrapolate_terminus_cells:
        lambda_field_description += (
            "; the grounding-line/grounded-marine-terminus band and "
            "every non-grounded cell (maskTerminusExtrapolated) were "
            "discarded and creep-fill extrapolated "
            f"({args.creep_fill_method}) from the grounded interior "
            "-- see --extrapolate-terminus-cells/"
            "creep_fill_extrapolate()"
        )

    out[args.lambda_field] = xr.DataArray(
        Lambda,
        dims=(ncell_dim,),
        attrs={
            "long_name": "Albany regularized-Coulomb bed roughness Lambda",
            "description": lambda_field_description,
        },
    )

    if args.diagnostics:
        if args.effective_pressure_type == "downs-johnson":
            effective_pressure_long_name = "Downs-Johnson effective pressure"
        elif args.effective_pressure_type == "ocean-connection":
            effective_pressure_long_name = (
                "Ocean-connection effective pressure (Hydrostatic "
                "Computed At Nodes, Use Pressurized Bed Above Sea "
                "Level: false)"
            )
        else:
            effective_pressure_long_name = (
                "Effective pressure (near-ocean/inland-transition "
                "parameterization; see effective_pressure4())"
            )

        out[args.effective_pressure_field] = xr.DataArray(
            N,
            dims=(ncell_dim,),
            attrs={
                "long_name": effective_pressure_long_name,
                "units": "Pa",
                "note": (
                    "N is a smooth function of bed elevation only "
                    "(not the actual flotation criterion, which also "
                    "depends on thickness), so some floating cells "
                    "near the grounding line (shallow bed) can have "
                    "nonzero N/floatationFraction != 1; this is "
                    "expected Albany behavior, not a bug."
                ),
            },
        )

        if args.source_effective_pressure_type == "constant":
            source_effective_pressure_long_name = (
                "Source-law effective pressure N_source (spatially "
                "uniform constant, kPa-equivalent scale)"
            )
            source_effective_pressure_units = "1 (kPa-equivalent)"
        elif args.source_effective_pressure_type == "field":
            source_effective_pressure_long_name = (
                f"Source-law effective pressure N_source (copied "
                f"from input field {args.source_effective_pressure_field}, "
                "converted to kPa-equivalent scale)"
            )
            source_effective_pressure_units = "1 (kPa-equivalent)"
        elif args.source_effective_pressure_type == "downs-johnson":
            source_effective_pressure_long_name = (
                "Source-law effective pressure N_source (Downs-"
                "Johnson parameterization, converted to "
                "kPa-equivalent scale)"
            )
            source_effective_pressure_units = "1 (kPa-equivalent)"
        elif args.source_effective_pressure_type == "ocean-connection":
            source_effective_pressure_long_name = (
                "Source-law effective pressure N_source "
                "(ocean-connection parameterization, converted to "
                "kPa-equivalent scale)"
            )
            source_effective_pressure_units = "1 (kPa-equivalent)"
        else:
            source_effective_pressure_long_name = (
                "Source-law effective pressure N_source (near-ocean/"
                "inland-transition parameterization, converted to "
                "kPa-equivalent scale)"
            )
            source_effective_pressure_units = "1 (kPa-equivalent)"

        out["effectivePressureSource"] = xr.DataArray(
            N_source_kpa,
            dims=(ncell_dim,),
            attrs={
                "long_name": source_effective_pressure_long_name,
                "units": source_effective_pressure_units,
                "description": (
                    "Effective pressure implied by the input source "
                    "friction law (Tau_b_source = mu * N_source * "
                    "speed^qW), independent of N/N_albany above -- "
                    "see --source-effective-pressure-type."
                ),
            },
        )

        # Floatation fraction: Pw / Pice, where Pw is basal water
        # pressure (Pice - N) and Pice is the ice overburden pressure
        # (rho_i * g * H). 0 at ice-free cells (Pice == 0).
        Pice = args.rho_ice * args.gravity * H
        Pw = Pice - N
        floatation_fraction = np.where(Pice > 0.0, Pw / np.where(Pice > 0.0, Pice, 1.0), 0.0)

        out[args.floatation_fraction_field] = xr.DataArray(
            floatation_fraction,
            dims=(ncell_dim,),
            attrs={
                "long_name": "Floatation fraction (Pw / Pice)",
                "description": (
                    "Pw = Pice - N (basal water pressure), Pice = "
                    "rho_i * gravity * thickness (ice overburden "
                    "pressure); 0 at ice-free cells."
                ),
                "units": "1",
            },
        )

        # Shreve hydraulic potential: phi = rho_freshwater * g * bed +
        # Pw. Uses freshwater density for the bed-elevation term
        # (subglacial water is fresh), distinct from the ocean/
        # seawater density (rho_water) used elsewhere in this script
        # for ocean-connectivity terms.
        hydropotential = args.rho_freshwater * args.gravity * bed + Pw

        out[args.hydropotential_field] = xr.DataArray(
            hydropotential,
            dims=(ncell_dim,),
            attrs={
                "long_name": "Shreve hydraulic potential",
                "description": (
                    "phi = rho_freshwater * gravity * bedTopography "
                    "+ Pw, with Pw = Pice - N"
                ),
                "units": "Pa",
            },
        )

        out["maskGrounded"] = xr.DataArray(
            grounded.astype(np.int8),
            dims=(ncell_dim,),
            attrs={
                "long_name": "Mask of grounded ice (Albany grounded-ice test)",
                "description": (
                    "1 where rho_i*H + rho_w*bed > 0 and H > 0, else 0"
                ),
            },
        )

        if args.method == "stress-match-fit":
            mask_fast_flowing_long_name = (
                "Mask of grounded, fast-flowing cells assumed to be in "
                "the full-Coulomb regime (used to fit C)"
            )
            mask_fast_flowing_description = (
                "1 where maskGrounded and speed (from velocity-x/y-"
                "field, last nVertInterfaces level) > critical "
                "velocity, else 0. bedRoughnessRC is forced to "
                f"{args.lambda_reference_value:g} at these cells."
            )
        else:
            mask_fast_flowing_long_name = (
                "Mask of grounded cells above the transition velocity "
                "u0 (diagnostic only; --method=transition-velocity "
                "does not use a fast/slow-region split)"
            )
            mask_fast_flowing_description = (
                "1 where maskGrounded and speed (from velocity-x/y-"
                "field, last nVertInterfaces level) > transition "
                "velocity u0 (--transition-velocity), else 0."
            )

        out["maskFastFlowing"] = xr.DataArray(
            fast_flowing.astype(np.int8),
            dims=(ncell_dim,),
            attrs={
                "long_name": mask_fast_flowing_long_name,
                "description": mask_fast_flowing_description,
            },
        )

        if args.method == "stress-match-fit":
            mask_valid_long_name = (
                "Mask of cells where bedRoughnessRC (Lambda) was "
                "solved exactly, rather than set to the full-Coulomb "
                f"reference value ({args.lambda_reference_value:g})"
            )
            mask_valid_description = (
                "1 where the cell is grounded, not fast-flowing, and "
                "the source-law basal shear stress at the cell's "
                "current sliding speed is strictly below the Coulomb "
                "limit C*N (a valid, non-negative Lambda solution "
                "exists); 0 otherwise (includes maskFastFlowing "
                "cells and any slow-flowing grounded cell where the "
                "solve is ill-defined)."
            )
        else:
            mask_valid_long_name = (
                "Mask of cells where bedRoughnessRC (Lambda) and "
                "muFriction (C) were computed exactly, rather than "
                f"set to their reference values "
                f"({args.lambda_reference_value:g}/"
                f"{args.mu_reference_value:g})"
            )
            mask_valid_description = (
                "1 where the cell is grounded and has a well-defined "
                "current sliding speed, N, and input mu (a valid "
                "closed-form solve exists -- see "
                "solve_transition_velocity()); 0 otherwise. Unlike "
                "--method=stress-match-fit, this is not restricted by "
                "maskFastFlowing -- every grounded cell with valid "
                "inputs gets an exact solve regardless of speed."
            )

        out["maskValidBedRoughnessRC"] = xr.DataArray(
            lambda_mask.astype(np.int8),
            dims=(ncell_dim,),
            attrs={
                "long_name": mask_valid_long_name,
                "description": mask_valid_description,
            },
        )

        if args.extrapolate_terminus_cells:
            out["maskTerminusExtrapolated"] = xr.DataArray(
                fill_mask.astype(np.int8),
                dims=(ncell_dim,),
                attrs={
                    "long_name": (
                        "Mask of the whole-domain region whose "
                        "computed bedRoughnessRC (and, for "
                        "--method=transition-velocity, muFriction) "
                        "values were discarded and replaced by "
                        "creep-fill extrapolation from the grounded "
                        "interior"
                    ),
                    "description": (
                        "1 where the cell is either a grounding-line/"
                        "grounded-marine-terminus cell (grounded and "
                        "adjacent to floating ice, or to ice-free, "
                        "bed-below-sea-level ocean) or any non-"
                        "grounded cell (floating ice, ice-free ocean, "
                        "ice-free land), else 0 (the grounded "
                        "interior, used as the extrapolation source). "
                        "See --extrapolate-terminus-cells/"
                        "--creep-fill-method/creep_fill_extrapolate()."
                    ),
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

    if args.diagnostics:
        out[args.flow_rate_field] = xr.DataArray(
            A,
            dims=(ncell_dim,),
            attrs={
                "long_name": flow_rate_long_name,
                "units": "Pa-3 s-1",
            },
        )

    # Save conversion information globally.
    out.attrs["regularizedCoulomb_method"] = args.method
    if args.method == "stress-match-fit":
        out.attrs["regularizedCoulomb_C"] = float(C)
        out.attrs["regularizedCoulomb_criticalVelocity"] = (
            float(args.critical_velocity)
        )
    else:
        out.attrs["regularizedCoulomb_transitionVelocity"] = (
            float(args.transition_velocity)
        )
        out.attrs["regularizedCoulomb_muReferenceValue"] = (
            float(args.mu_reference_value)
        )
    out.attrs["regularizedCoulomb_q"] = float(RC_POWER_EXPONENT)
    out.attrs["weertman_q"] = float(args.weertman_q)
    out.attrs["regularizedCoulomb_GlenN"] = float(args.glen_n)
    out.attrs["regularizedCoulomb_flowRateType"] = albany_flow_rate_type
    out.attrs["regularizedCoulomb_effectivePressureType"] = (
        args.effective_pressure_type
    )
    if args.effective_pressure_type != "ocean-connection":
        out.attrs["regularizedCoulomb_minFractionOverburden"] = (
            float(args.min_fraction_overburden)
        )
        out.attrs["regularizedCoulomb_pressureLengthScale"] = (
            float(args.pressure_length_scale)
        )
    if args.effective_pressure_type == "transition":
        out.attrs["regularizedCoulomb_transitionHOcean"] = (
            float(args.transition_h_ocean)
        )
    out.attrs["regularizedCoulomb_sourceEffectivePressureType"] = (
        args.source_effective_pressure_type
    )
    if args.source_effective_pressure_type == "constant":
        out.attrs["regularizedCoulomb_sourceEffectivePressure"] = (
            float(args.source_effective_pressure)
        )
    elif args.source_effective_pressure_type == "field":
        out.attrs["regularizedCoulomb_sourceEffectivePressureField"] = (
            args.source_effective_pressure_field
        )
    if args.source_effective_pressure_type in ("downs-johnson", "transition"):
        out.attrs["regularizedCoulomb_sourceMinFractionOverburden"] = (
            float(args.source_min_fraction_overburden)
        )
        out.attrs["regularizedCoulomb_sourcePressureLengthScale"] = (
            float(args.source_pressure_length_scale)
        )
    if args.source_effective_pressure_type == "transition":
        out.attrs["regularizedCoulomb_sourceTransitionHOcean"] = (
            float(args.source_transition_h_ocean)
        )

    # xarray cannot safely overwrite an open source file, so use temp.
    tmp = args.output + ".tmp"
    out.to_netcdf(tmp)
    out.close()

    shutil.move(tmp, args.output)

    print(f"Wrote converted IC: {args.output}")

    if args.plot_transects:
        x_cell = np.asarray(ds["xCell"].values, dtype=np.float64)
        y_cell = np.asarray(ds["yCell"].values, dtype=np.float64)

        transect_fields = {
            "N": (N, "Pa", effective_pressure_long_name),
            "floatation fraction": (
                floatation_fraction, "1", "Floatation fraction (Pw / Pice)"
            ),
            "hydropotential": (
                hydropotential, "Pa", "Shreve hydraulic potential"
            ),
        }
        # "implied C" is only meaningful (and only computed) for
        # --method=stress-match-fit -- --method=transition-velocity
        # already writes an exactly-matching, spatially-varying C
        # directly to --mu-field.
        if args.method == "stress-match-fit":
            transect_fields["implied C"] = (
                local_C, "1",
                "Implied local C (source Tau_b / N), Coulomb "
                "C-fit region only",
            )

        plot_transects(
            transect_names=args.transect_names,
            transects_dir=args.transects_dir,
            plot_dir=args.plot_dir,
            x_cell=x_cell,
            y_cell=y_cell,
            fields=transect_fields,
            thickness=H,
            bed=bed,
            rho_i=args.rho_ice,
            rho_w=args.rho_water,
            min_fraction_overburden=args.min_fraction_overburden,
            fit_mask=fit_mask if args.method == "stress-match-fit" else None,
            fitted_c=C if args.method == "stress-match-fit" else None,
        )

        # Antarctic-wide maps of the same diagnostic fields plus
        # bedRoughnessRC (alongside the original Weertman muFriction
        # field for comparison, both log scale), using mosaic
        # (e3sm-unified environment). Always produced alongside the
        # transect plots -- no separate CLI flag -- since both are
        # diagnostic/evaluation outputs triggered by --plot-transects.
        map_fields = {
            args.lambda_field: [
                {
                    "values": Lambda, "units": "m",
                    "title": (
                        "Albany regularized-Coulomb bed roughness Lambda"
                    ),
                    "log": True, "cmap": "turbo",
                },
            ],
            args.mu_field: [
                {
                    "values": mu,
                    "units": (
                        f"kPa (m yr-1)^-{args.weertman_q:g}"
                        if args.source_effective_pressure_type == "constant"
                        else f"(m yr-1)^-{args.weertman_q:g}"
                    ),
                    "title": "Original source-law muFriction (input)",
                    "log": True, "cmap": "turbo_r",
                },
                {
                    "values": mu_field_values,
                    "units": "1",
                    "title": (
                        f"Output {args.mu_field} (Regularized Coulomb C)"
                    ),
                    "log": True, "cmap": "turbo",
                },
            ],
            args.effective_pressure_field: [
                {
                    "values": N, "units": "Pa",
                    "title": effective_pressure_long_name,
                    # N is only physically meaningful for grounded
                    # ice (see the floating-cell-nonzero-N note on
                    # this field in the output NetCDF).
                    "mask": grounded,
                }
            ],
            args.floatation_fraction_field: [
                {
                    "values": floatation_fraction, "units": "1",
                    "title": "Floatation fraction (Pw / Pice)",
                    # A few outlier cells (very thin/no ice) can
                    # produce huge ratios; clip the colormap to the
                    # physically meaningful [0, 1] range, and mask
                    # out non-grounded cells entirely (this field is
                    # only meaningful for grounded ice).
                    "vmax": 1.0, "mask": grounded,
                }
            ],
            args.hydropotential_field: [
                {
                    "values": hydropotential, "units": "Pa",
                    "title": "Shreve hydraulic potential",
                    "mask": grounded,
                }
            ],
            args.flow_rate_field: [
                {
                    "values": A, "units": "Pa-3 s-1",
                    "title": flow_rate_long_name,
                }
            ],
            "maskGrounded": [
                {
                    "values": grounded.astype(np.int8), "units": "1",
                    "title": "Mask of grounded ice",
                }
            ],
            "maskFastFlowing": [
                {
                    "values": fast_flowing.astype(np.int8), "units": "1",
                    "title": "Mask of grounded, fast-flowing cells",
                }
            ],
            "maskValidBedRoughnessRC": [
                {
                    "values": lambda_mask.astype(np.int8), "units": "1",
                    "title": (
                        "Mask of cells where bedRoughnessRC was solved "
                        "exactly"
                    ),
                }
            ],
            "regimeRatio": [
                {
                    "values": regime_ratio, "units": "1",
                    "title": (
                        "u / (Lambda * A * N^n): Coulomb (>1, red) vs. "
                        "power-law (<1, blue) regime"
                    ),
                    "log": True, "cmap": "RdBu_r",
                    "vmin": 0.1, "vmax": 10.0,
                    "mask": grounded,
                }
            ],
        }

        # "impliedC" is only meaningful (and only computed) for
        # --method=stress-match-fit -- --method=transition-velocity
        # already writes an exactly-matching, spatially-varying C
        # directly to --mu-field (visible via the "muFriction (input)"
        # panel above, which shows the *original* input mu instead).
        if args.method == "stress-match-fit":
            map_fields["impliedC"] = [
                {
                    "values": local_C, "units": "1",
                    "title": (
                        "Implied C (source Tau_b / N) in the fast-"
                        f"flowing/full-Coulomb fit region; fitted "
                        f"scalar C = {C:.4g}"
                    ),
                    "log": True, "cmap": "turbo",
                    "mask": fit_mask,
                }
            ]

        map_transects = None
        if args.plot_transects_on_maps:
            import pyproj

            transformer = pyproj.Transformer.from_crs(
                "epsg:4326", "epsg:3031", always_xy=True
            )
            map_transects = []
            for name in args.transect_names:
                t_lon, t_lat = load_transect(name, args.transects_dir)
                t_x, t_y, _ = project_transect(t_lon, t_lat, transformer)
                map_transects.append((name, t_x, t_y))

        plot_maps(
            mesh_ds=ds,
            fields=map_fields,
            plot_dir=args.plot_dir,
            transects=map_transects,
        )

    if args.flow_rate_type == "constant":
        flow_rate_yaml_lines = (
            f"        Flow Rate Type: Constant\n"
            f"        Flow Rate: {args.flow_rate:.16e}\n"
        )
    else:
        flow_rate_yaml_lines = (
            f"        Flow Rate Type: Temperature Based\n"
        )

    if args.effective_pressure_type == "downs-johnson":
        effective_pressure_yaml_lines = (
            "        Effective Pressure Type: Hydrostatic Computed At "
            "Nodes\n"
            "        Use Pressurized Bed Above Sea Level: true\n"
            "        Minimum Fraction Overburden Pressure: "
            f"{args.min_fraction_overburden:.16e}\n"
            "        Length Scale Factor: "
            f"{args.pressure_length_scale / 1000.0:.16e}"
        )
    elif args.effective_pressure_type == "ocean-connection":
        effective_pressure_yaml_lines = (
            "        Effective Pressure Type: Hydrostatic Computed At "
            "Nodes\n"
            "        Use Pressurized Bed Above Sea Level: false"
        )
    else:
        effective_pressure_yaml_lines = (
            "        Effective Pressure Type: TRANSITION OPTION TO BE ADDED"
        )

    print()
    print("Suggested Albany YAML section")
    print("-----------------------------")
    print(
    f"""
      Basal Friction Coefficient:
        Type: Regularized Coulomb
        Mu Type: Field
        Power Exponent: {RC_POWER_EXPONENT:.16e}
{flow_rate_yaml_lines}
        Bed Roughness Type: Field
{effective_pressure_yaml_lines}
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
    if args.effective_pressure_type == "transition":
        print(
            "NOTE: The \"transition\" effective-pressure "
            "parameterization (--effective-pressure-type=transition) "
            "is not yet implemented in Albany, so no valid "
            "\"Effective Pressure Type\" YAML setting exists for it "
            "yet -- the placeholder above must be replaced once "
            "Albany supports reading a precomputed N field directly "
            f"(e.g. from the \"{args.effective_pressure_field}\" field "
            f"written to {args.output}) for this parameterization."
        )



if __name__ == "__main__":
    main()
