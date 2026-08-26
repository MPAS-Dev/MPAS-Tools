#!/usr/bin/env python3
"""
Tangential "what-if" diagnostic for the Weertman -> Regularized
Coulomb (RC) conversion carried out by friction_law_conversion.py.

For a user-supplied list of candidate scalar Coulomb Friction
Coefficients C, plot the effective pressure N that *pure* Coulomb
sliding (no regularization: Tau_b = C * N, i.e. the RC law's
fully-plastic limit) would imply at each cell's actual current basal
shear stress:

    Tau_b_weertman = mu * speed^qW      (MALI's Weertman law)
    N(C)           = Tau_b_weertman / C

This is plotted along the same glacier transects used by
friction_law_conversion.py's --plot-transects, restricted to the same
"fast-flowing" region (grounded cells with speed > critical velocity)
that script uses to fit its own single scalar C -- i.e. exactly the
region where the pure-Coulomb assumption is invoked -- so the curves
can be visually compared against a range of candidate C values,
including (optionally) the one friction_law_conversion.py itself
would fit.

This script is intentionally lightweight: it reuses
friction_law_conversion.py's transect-loading/projection/plotting
machinery directly (via import) rather than duplicating it, and
duplicates only the handful of lines of physics (mu * speed^qW) needed
to get from the MALI input file to Tau_b_weertman.

Example
-------
    python3 coulomb_N_sensitivity.py relaxed_10yrs_4km.nc \\
        --uc 100 --weertman-q 0.2 \\
        --transects-dir /path/to/geometric_data/landice/transect \\
        --c-values 0.005 0.01 0.0162 0.02 0.05 \\
        --plot-dir coulomb_N_sensitivity_plots
"""

import argparse

import numpy as np
import xarray as xr

from friction_law_conversion import (
    SECONDS_PER_YEAR,
    plot_transects,
)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", help="Input MALI initial-condition NetCDF file")

    parser.add_argument(
        "--critical-velocity", "--uc",
        type=float,
        required=True,
        help=(
            "Critical velocity u_c [m/yr]: cells with a current "
            "sliding speed above this, within grounded ice, define "
            "the fast-flowing region the pure-Coulomb N(C) curves "
            "are restricted to (same convention as "
            "friction_law_conversion.py)."
        ),
    )
    parser.add_argument(
        "--weertman-q", "--q",
        dest="weertman_q",
        type=float,
        default=0.2,
        help="Input Weertman/Power-Law sliding exponent qW (default: 0.2)",
    )
    parser.add_argument(
        "--c-values",
        type=float,
        nargs="+",
        required=True,
        help=(
            "List of candidate scalar Coulomb Friction Coefficients C "
            "to evaluate, e.g. --c-values 0.005 0.01 0.02."
        ),
    )

    parser.add_argument("--rho-ice", type=float, default=910.0)
    parser.add_argument("--rho-water", type=float, default=1028.0)
    parser.add_argument(
        "--rho-freshwater", type=float, default=1000.0,
        help=(
            "Freshwater density [kg m^-3] (default: 1000.0), used "
            "for the bed-elevation term of the hydropotential panel "
            "(same convention as friction_law_conversion.py)."
        ),
    )
    parser.add_argument("--gravity", type=float, default=9.80616)

    parser.add_argument(
        "--mu-field", default="muFriction",
        help="Weertman friction field (default: muFriction)",
    )
    parser.add_argument(
        "--thickness-field", default="thickness",
        help="Ice thickness field (default: thickness)",
    )
    parser.add_argument(
        "--bed-field", default="bedTopography",
        help="Bed elevation field (default: bedTopography)",
    )
    parser.add_argument(
        "--velocity-x-field", default="uReconstructX",
        help="MALI x-velocity field [m s^-1] (default: uReconstructX)",
    )
    parser.add_argument(
        "--velocity-y-field", default="uReconstructY",
        help="MALI y-velocity field [m s^-1] (default: uReconstructY)",
    )
    parser.add_argument(
        "--time-index", type=int, default=0,
        help="Time index for Time-dependent IC fields (default: 0)",
    )

    parser.add_argument(
        "--transects-dir",
        required=True,
        help=(
            "Path to a geometric_features landice/transect directory "
            "(see friction_law_conversion.py --transects-dir)."
        ),
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
            "--transects-dir)."
        ),
    )
    parser.add_argument(
        "--plot-dir",
        default="coulomb_N_sensitivity_plots",
        help=(
            "Directory to write output PNGs to (default: "
            "coulomb_N_sensitivity_plots, created if needed)."
        ),
    )

    args = parser.parse_args()

    ds = xr.open_dataset(args.input, mask_and_scale=False)

    def cell_field(name):
        da = ds[name]
        if "Time" in da.dims:
            da = da.isel(Time=args.time_index)
        vert_dims = [d for d in da.dims if d.lower().startswith("nvert")]
        if vert_dims:
            da = da.isel({vert_dims[0]: -1})
        return np.asarray(da.values).squeeze().astype(np.float64)

    mu = cell_field(args.mu_field)
    H = cell_field(args.thickness_field)
    bed = cell_field(args.bed_field)
    uX = cell_field(args.velocity_x_field)
    uY = cell_field(args.velocity_y_field)
    x_cell = np.asarray(ds["xCell"].values, dtype=np.float64)
    y_cell = np.asarray(ds["yCell"].values, dtype=np.float64)

    # Basal sliding speed [m/yr], matching Albany's internal
    # convention (see friction_law_conversion.py).
    speed = np.sqrt(uX ** 2 + uY ** 2) * SECONDS_PER_YEAR

    # MALI's Weertman sliding law has no effective-pressure term:
    # Tau_b = mu * speed^qW (same as fit_coulomb_C_fast_region() in
    # friction_law_conversion.py).
    tau_b_weertman = mu * speed ** args.weertman_q

    # Same grounded-ice/fast-flowing tests friction_law_conversion.py
    # uses to define the region assumed to already be in the
    # fully-plastic Coulomb regime.
    grounded = (H > 0.0) & (args.rho_ice * H + args.rho_water * bed > 0.0)
    fast_flowing = grounded & (speed > args.critical_velocity)

    # Pure-Coulomb inversion: Tau_b = C * N  =>  N(C) = Tau_b / C, in
    # physical Pa (no Albany-internal kPa rescaling needed here since
    # this script works entirely in physical units). For each
    # candidate C, also derive the floatation fraction and
    # hydropotential that N(C) would imply, using the same formulas
    # as friction_law_conversion.py (Pice = rho_i * g * H,
    # Pw = Pice - N, floatation_fraction = Pw / Pice, hydropotential =
    # rho_freshwater * g * bed + Pw), all restricted to the
    # fast-flowing region.
    Pice = args.rho_ice * args.gravity * H

    n_of_c = {}
    floatation_fraction_of_c = {}
    hydropotential_of_c = {}
    for c in args.c_values:
        label = f"C={c:g}"

        n_c = np.full_like(tau_b_weertman, np.nan)
        n_c[fast_flowing] = tau_b_weertman[fast_flowing] / c
        n_of_c[label] = n_c

        Pw_c = Pice - n_c
        floatation_fraction_of_c[label] = np.where(
            fast_flowing & (Pice > 0.0), Pw_c / Pice, np.nan
        )
        hydropotential_of_c[label] = np.where(
            fast_flowing, args.rho_freshwater * args.gravity * bed + Pw_c,
            np.nan,
        )

    plot_transects(
        transect_names=args.transect_names,
        transects_dir=args.transects_dir,
        plot_dir=args.plot_dir,
        x_cell=x_cell,
        y_cell=y_cell,
        fields={
            "N(C)": (
                n_of_c, "Pa",
                "Pure-Coulomb-implied N for candidate C values "
                "(fast-flowing region only)",
            ),
            "floatation fraction(C)": (
                floatation_fraction_of_c, "1",
                "Pure-Coulomb-implied floatation fraction (Pw / Pice) "
                "for candidate C values (fast-flowing region only)",
                (0.0, 1.0),
            ),
            "hydropotential(C)": (
                hydropotential_of_c, "Pa",
                "Pure-Coulomb-implied Shreve hydraulic potential for "
                "candidate C values (fast-flowing region only)",
            ),
        },
        thickness=H,
        bed=bed,
        rho_i=args.rho_ice,
        rho_w=args.rho_water,
        fit_mask=fast_flowing,
    )


if __name__ == "__main__":
    main()
