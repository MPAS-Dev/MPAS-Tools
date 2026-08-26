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
    mu       : kPa * (yr/m)^qW (per a correction to MPAS-Tools'
               Registry.xml; NOT Pa * (yr/m)^qW), so that
               mu * speed^qW comes out already in the same kPa scale
               as N_albany (see ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT
               and fit_coulomb_C_fast_region())
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


def fit_coulomb_C_fast_region(mu, N, area, speed, q, mask):
    """
    Fit a single scalar Regularized Coulomb "Coulomb Friction
    Coefficient" C by assuming that, in the fast-flowing region
    identified by `mask` (typically grounded cells with a current
    sliding speed above some critical velocity), the ice is already in
    the fully-plastic Coulomb regime of the RC law, i.e.

        Tau_b_RC = C * N

    C is chosen to be the area-weighted mean, over the fast-flowing
    region, of the per-cell "local C" implied by that relation against
    the Weertman law's basal shear stress at the cell's actual current
    sliding speed,

        Tau_b_Weertman = mu * speed^qW

    (no effective-pressure term -- MALI's Weertman sliding law has
    none; muFriction was calibrated against this convention).

        local_C = Tau_b_Weertman / N

        C = sum(area * local_C) / sum(area)

    A straight area-weighted mean is used (rather than an N^2-weighted
    least-squares fit through the origin) so that a small number of
    high-N, high-leverage cells cannot dominate the result.

    `q` here is the input Weertman/Power-Law exponent (qW), not the
    Regularized Coulomb exponent.

    IMPORTANT: `N` must already be expressed in the same units Albany
    will actually use at runtime for its own internal effective
    pressure (numerically equal to physical Pa / 1000, i.e. "kPa";
    see ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT), not raw SI Pascals.
    Passing raw-Pa N here would make C inconsistent with how Albany
    multiplies C against its own internal N at runtime. This is
    dimensionally consistent with `mu` (MALI's `muFriction`), whose
    correct physical units -- per a correction to MPAS-Tools'
    Registry.xml -- are kPa * (yr/m)^qW (not Pa * (yr/m)^qW as one
    might otherwise assume): Tau_b_Weertman = mu * speed^qW therefore
    comes out in kPa already, matching `N`'s kPa scale here, with no
    separate unit conversion needed in this function.
    """
    valid = (
        mask
        & np.isfinite(mu)
        & np.isfinite(N)
        & np.isfinite(area)
        & np.isfinite(speed)
        & (N > 0.0)
        & (area > 0.0)
        & (speed > 0.0)
    )

    if not np.any(valid):
        raise ValueError(
            "No valid fast-flowing cells (speed > critical velocity) "
            "available for C calculation."
        )

    tau_b_weertman = mu[valid] * speed[valid] ** q
    local_C = tau_b_weertman / N[valid]

    C = np.sum(area[valid] * local_C) / np.sum(area[valid])

    return C, valid


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
                out-of-range outliers); ignored if `log` is True.
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
                    vmin=np.nanmin(plot_values),
                    vmax=np.nanmax(plot_values),
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

    parser.add_argument(
        "--lambda-reference-value",
        type=float,
        default=0.0,
        help=(
            "Value assigned to bedRoughnessRC (Lambda) at cells "
            "assumed to be in the full-Coulomb regime: all "
            "fast-flowing cells (speed > critical velocity, the same "
            "region used to fit C), plus any other grounded cell "
            "where the exact per-cell Lambda solve is ill-defined "
            "(Weertman Tau_b already meets or exceeds the Coulomb "
            "limit C*N). Lambda -> 0 exactly reproduces the "
            "full-Coulomb limit, so 0.0 (default) is physically "
            "correct; a small positive reference value can be used "
            "instead if a strictly-zero bed roughness is undesirable "
            "for other reasons (default: 0.0)."
        )
    )

    # Effective pressure N
    parser.add_argument(
        "--effective-pressure-type",
        choices=["downs-johnson", "transition"],
        default="downs-johnson",
        help=(
            "How to compute the effective pressure N (default: "
            "downs-johnson). \"downs-johnson\" reproduces Albany's "
            "own internal \"Hydrostatic At Nodes\" Effective Pressure "
            "Type formula exactly (see --min-fraction-overburden/"
            "--pressure-length-scale); Albany recomputes N itself at "
            "runtime from thickness/bed, so no N field needs to be "
            "supplied. \"transition\" computes N here using a "
            "near-ocean/inland-transition parameterization (see "
            "--min-fraction-overburden/--pressure-length-scale/"
            "--transition-h-ocean) and writes it to the output for "
            "reference; NOTE: Albany does not yet have a way to "
            "consume this precomputed N directly for the Regularized "
            "Coulomb law, so \"transition\" cannot currently be used "
            "to actually run Albany (offline evaluation only, pending "
            "upstream Albany support)."
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
            "Weertman friction field (default: muFriction). Its "
            "correct physical units -- per a correction to "
            "MPAS-Tools' Registry.xml -- are kPa * (yr/m)^qW (not "
            "Pa * (yr/m)^qW); this script's use of mu is dimensionally "
            "consistent with that (see fit_coulomb_C_fast_region())."
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
            "production runs where only bedRoughnessRC is needed and "
            "extra fields are unwanted clutter."
        )
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

    if args.min_fraction_overburden is None or args.pressure_length_scale is None:
        parser.error(
            "--min-fraction-overburden and --pressure-length-scale "
            "are required (used by both "
            "--effective-pressure-type=downs-johnson and "
            "--effective-pressure-type=transition)"
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
    # Optimal C: fit against the fast-flowing region only, assuming
    # it is already in the fully-plastic Coulomb regime of the RC law
    # (Tau_b = C * N).
    # -------------------------------------------------------------
    fast_flowing = grounded & (speed > args.critical_velocity)

    C, fit_mask = fit_coulomb_C_fast_region(
        mu=mu,
        N=N_albany,
        area=area,
        speed=speed,
        q=args.weertman_q,
        mask=fast_flowing,
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
    # This is the same Weertman shear-stress convention used by
    # fit_coulomb_C_fast_region above (Tau_b_W = mu * u^qW, no N term,
    # matched in the fast-flowing region against the RC Coulomb limit
    # C*N).
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
    # (non-negative) solution for Lambda. Physically, Lambda -> 0 is
    # *exactly* the fully-plastic Coulomb regime (Tau_b_RC saturates
    # at its maximum achievable value, C*N, independent of speed), so
    # setting Lambda = args.lambda_reference_value at these cells is
    # not an arbitrary filler value -- it is the correct behavior for
    # cells that the fast-flowing/full-Coulomb assumption is designed
    # to describe in the first place.
    #
    # Fast-flowing cells (the same region used to fit C, i.e.
    # speed > critical_velocity) are *always* assumed to be in this
    # full-Coulomb regime and are therefore always forced to
    # Lambda = args.lambda_reference_value, regardless of what the
    # per-cell algebraic solve above would otherwise give -- the exact
    # per-cell solve is not attempted there at all, since matching the
    # Weertman law exactly at high speed is not the goal (the fast
    # region is assumed C-limited by construction).
    # -------------------------------------------------------------
    Lambda = np.full_like(N, args.lambda_reference_value)

    # tau_b_weertman/local_C are computed over all grounded cells with
    # a well-defined speed and mu (independent of the fast/slow split)
    # so that diagnostics (local_C) remain meaningful for the
    # fast-flowing fit region even though the exact Lambda solve below
    # is only attempted for the slow-flowing cells.
    speed_defined = (
        grounded
        & np.isfinite(N) & (N > 0.0)
        & np.isfinite(mu) & (mu > 0.0)
        & np.isfinite(speed) & (speed > 0.0)
    )

    tau_b_weertman = np.full_like(N, np.nan)
    tau_b_weertman[speed_defined] = (
        mu[speed_defined] * speed[speed_defined] ** args.weertman_q
    )

    valid_speed = speed_defined & ~fast_flowing

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
            f"NOTE: {n_unreachable} slow-flowing grounded cells have a "
            "Weertman basal shear stress (mu*u^qW) at the current "
            "sliding speed that meets or exceeds the Coulomb limit "
            f"C*N; Lambda set to {args.lambda_reference_value:g} "
            "(maximal Coulomb sliding) at these cells."
        )
    print(
        f"Fast-flowing cells forced to Lambda = "
        f"{args.lambda_reference_value:g} (full-Coulomb assumption) : "
        f"{int(np.count_nonzero(fast_flowing))}"
    )

    # Floating/ice-free/invalid-speed/unreachable-stress cells are
    # deliberately left at zero (see Lambda initialization above).

    # -------------------------------------------------------------
    # Diagnostics
    # -------------------------------------------------------------
    # "Local" implied C: the per-cell ratio of the actual Weertman
    # shear stress (at the cell's current sliding speed) to N, i.e.
    # what C would have to be for that cell alone to be exactly in
    # the full-Coulomb regime. Its spread within the fast-flowing fit
    # region gives a sense of how well a single scalar C fits that
    # region.
    local_C = np.full_like(N, np.nan)
    local_C[speed_defined] = (
        tau_b_weertman[speed_defined] / N_albany[speed_defined]
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
                "Fast-flowing cells (speed > critical velocity) and "
                "any other grounded cell where the solve below is "
                "ill-defined are assumed to be in the full-Coulomb "
                f"regime and set to {args.lambda_reference_value:g} "
                "(see maskFastFlowing/maskValidBedRoughnessRC). "
                "Elsewhere, Lambda is solved exactly so that the "
                "Regularized Coulomb law reproduces the Weertman "
                "law's basal shear stress (mu*u^qW, no effective-"
                "pressure term) at the cell's actual current sliding "
                "speed u (from velocity-x/y-field, last "
                "nVertInterfaces level): Lambda = u * "
                "[(C*N/(mu*u^qW))^(1/qR) - 1] / (SECONDS_PER_YEAR * A "
                "* N^n), with u in m/yr, A in Pa^-3 s^-1, N in Pa, "
                "matching Albany's internal secsInYr scaling in "
                "LandIce_BasalFrictionCoefficient_Def.hpp"
            ),
        },
    )

    if args.diagnostics:
        if args.effective_pressure_type == "downs-johnson":
            effective_pressure_long_name = "Downs-Johnson effective pressure"
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

        out["maskFastFlowing"] = xr.DataArray(
            fast_flowing.astype(np.int8),
            dims=(ncell_dim,),
            attrs={
                "long_name": (
                    "Mask of grounded, fast-flowing cells assumed to be in "
                    "the full-Coulomb regime (used to fit C)"
                ),
                "description": (
                    "1 where maskGrounded and speed (from velocity-x/y-"
                    "field, last nVertInterfaces level) > critical "
                    "velocity, else 0. bedRoughnessRC is forced to "
                    f"{args.lambda_reference_value:g} at these cells."
                ),
            },
        )

        out["maskValidBedRoughnessRC"] = xr.DataArray(
            lambda_mask.astype(np.int8),
            dims=(ncell_dim,),
            attrs={
                "long_name": (
                    "Mask of cells where bedRoughnessRC (Lambda) was "
                    "solved exactly, rather than set to the full-Coulomb "
                    f"reference value ({args.lambda_reference_value:g})"
                ),
                "description": (
                    "1 where the cell is grounded, not fast-flowing, and "
                    "the Weertman basal shear stress at the cell's "
                    "current sliding speed is strictly below the Coulomb "
                    "limit C*N (a valid, non-negative Lambda solution "
                    "exists); 0 otherwise (includes maskFastFlowing "
                    "cells and any slow-flowing grounded cell where the "
                    "solve is ill-defined)."
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
    out.attrs["regularizedCoulomb_C"] = float(C)
    out.attrs["regularizedCoulomb_criticalVelocity"] = (
        float(args.critical_velocity)
    )
    out.attrs["regularizedCoulomb_q"] = float(RC_POWER_EXPONENT)
    out.attrs["weertman_q"] = float(args.weertman_q)
    out.attrs["regularizedCoulomb_GlenN"] = float(args.glen_n)
    out.attrs["regularizedCoulomb_flowRateType"] = albany_flow_rate_type
    out.attrs["regularizedCoulomb_effectivePressureType"] = (
        args.effective_pressure_type
    )
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

    # xarray cannot safely overwrite an open source file, so use temp.
    tmp = args.output + ".tmp"
    out.to_netcdf(tmp)
    out.close()

    shutil.move(tmp, args.output)

    print(f"Wrote converted IC: {args.output}")

    if args.plot_transects:
        x_cell = np.asarray(ds["xCell"].values, dtype=np.float64)
        y_cell = np.asarray(ds["yCell"].values, dtype=np.float64)

        plot_transects(
            transect_names=args.transect_names,
            transects_dir=args.transects_dir,
            plot_dir=args.plot_dir,
            x_cell=x_cell,
            y_cell=y_cell,
            fields={
                "N": (N, "Pa", effective_pressure_long_name),
                "floatation fraction": (
                    floatation_fraction, "1", "Floatation fraction (Pw / Pice)"
                ),
                "hydropotential": (
                    hydropotential, "Pa", "Shreve hydraulic potential"
                ),
                "implied C": (
                    local_C, "1",
                    "Implied local C (Weertman Tau_b / N), Coulomb "
                    "C-fit region only",
                ),
            },
            thickness=H,
            bed=bed,
            rho_i=args.rho_ice,
            rho_w=args.rho_water,
            min_fraction_overburden=args.min_fraction_overburden,
            fit_mask=fit_mask,
            fitted_c=C,
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
                    "values": Lambda, "units": "Pa (m yr-1)^-1/3",
                    "title": (
                        "Albany regularized-Coulomb bed roughness Lambda"
                    ),
                    "log": True, "cmap": "turbo",
                },
                {
                    "values": mu,
                    "units": f"kPa (m yr-1)^-{args.weertman_q:g}",
                    "title": "Original Weertman muFriction (input)",
                    "log": True, "cmap": "turbo_r",
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
            "impliedC": [
                {
                    "values": local_C, "units": "1",
                    "title": (
                        "Implied C (Weertman Tau_b / N) in the fast-"
                        f"flowing/full-Coulomb fit region; fitted "
                        f"scalar C = {C:.4g}"
                    ),
                    "log": True, "cmap": "turbo",
                    "mask": fit_mask,
                }
            ],
        }

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
            "        Effective Pressure Type: Hydrostatic At Nodes\n"
            "        Use Pressurized Bed Above Sea Level: true\n"
            "        Minimum Fraction Overburden Pressure: "
            f"{args.min_fraction_overburden:.16e}\n"
            "        Length Scale Factor: "
            f"{args.pressure_length_scale / 1000.0:.16e}"
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
    LandIce BCs:
      Basal Friction Coefficient:
        Type: Regularized Coulomb
        Coulomb Friction Coefficient: {C:.16e}
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
