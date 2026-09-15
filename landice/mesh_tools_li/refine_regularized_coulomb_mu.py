#!/usr/bin/env python3
"""
Refine a Regularized Coulomb (RC) friction-law C (muFriction) field
via a fixed-point (Picard) update against Albany's own actual
simulated velocity, as an alternative to friction_law_conversion.py's
--fe-correct-mu.

Workflow
--------
1. Run friction_law_conversion.py once, as usual, on the original
   MALI initial condition to produce a Lambda/C conversion, e.g.:

       friction_law_conversion.py original_input.nc rc_pass1.nc ...

2. Configure Albany's YAML to read Lambda/C from rc_pass1.nc, run it,
   and obtain a MALI-compatible output file containing the *actual*
   resulting basal velocity (e.g. `albany_pass1_output.nc`, from
   however you already convert/translate Albany's Exodus output back
   to MPAS format).
3. Run this script:

       refine_regularized_coulomb_mu.py original_input.nc rc_pass1.nc \\
           albany_pass1_output.nc rc_pass2.nc [options]

   to produce a refined `rc_pass2.nc`, with C updated and Lambda left
   unchanged. All physics options (weertman-q, effective-pressure-
   type, flow-rate, etc.) are read automatically from the
   `rc_pass1.nc.config.json` sidecar file that friction_law_
   conversion.py writes next to its output, so they do not need to be
   (and should not be) re-specified by hand here -- this avoids a
   common, hard-to-spot source of error where a physics option given
   to this script silently differs from the one used for the
   original conversion. See --config for details, including how to
   override individual options when intentionally desired.
4. Re-run Albany with rc_pass2.nc in place of rc_pass1.nc, and repeat
   from step 2 as needed -- rc_pass2.nc becomes the new "previous
   conversion output" and its own Albany output becomes the new
   "model output" -- watching the printed velocity-mismatch summary
   each time to judge convergence.

Why this can do better than --fe-correct-mu
--------------------------------------------
--fe-correct-mu only *approximates* how Albany's FEM evaluates the
friction law, by reconstructing Albany's own quadrature-point
interpolation of N/C/Lambda/speed from MPAS mesh connectivity. This
script instead uses Albany's *actual* simulated velocity, so it
automatically reflects everything Albany really does -- the true FE
quadrature-point evaluation, mesh/element details, and the full
nonlinear force-balance coupling between neighboring cells -- not
just an approximation of the friction law evaluated in isolation.

The update
----------
Holding N, N_albany, and Lambda fixed at their original
friction_law_conversion.py values (only C is refined; see that
script's --fe-correct-mu module docstring for why C alone is
targeted for this kind of correction), the Regularized Coulomb
stress is exactly linear in C:

    Tau_b_RC(u) = C * k(u),   k(u) = N_albany * u^qR / (u + u_c)^qR

Given Albany's actual simulated velocity u_model (from step 2 above),
the C that makes the RC law reproduce, AT THAT SAME VELOCITY u_model,
the stress the *original* (trusted) source friction law would predict
there --

    Tau_b_target(u_model) = mu * N_source * u_model^qW

(mu, N_source, and qW all come from the *original* input file/CLI
options, exactly as in the very first friction_law_conversion.py run
that solved Lambda/C -- but evaluated at u_model, Albany's own
simulated speed, NOT at the original input file's velocity) -- is

    C_new = Tau_b_target(u_model) / k(u_model)

(this is exactly friction_law_conversion.py's own
solve_transition_velocity() C formula, re-solved with u_model in
place of the original conversion's speed). This is applied as an
(optionally relaxed) multiplicative Picard update:

    C_refined = C_old * (C_new / C_old) ** relaxation

with relaxation=1 (the default) giving C_new exactly, and
relaxation<1 damping the step -- useful if repeated iterations
oscillate rather than converge, since Albany's actual velocity
response to a change in C is not perfectly local (membrane stresses
couple neighboring cells, so a purely pointwise update is only an
approximation of the true, globally-coupled sensitivity).

IMPORTANT: Tau_b_target must be evaluated at u_model, not at the
original input file's velocity. k(u) is monotonically increasing in
u (saturating at N_albany as u -> infinity), so pinning Tau_b_target
at a fixed velocity while evaluating k(.) at a *different* u_model
is structurally backwards: whenever u_model is larger than that fixed
velocity, it forces C to *decrease* (less friction), making an
already-too-fast region even faster -- and the reverse in regions
where the RC law was already close, needlessly perturbing them. Only
evaluating both Tau_b_target and k(.) at the *same* velocity, u_model,
gives an update that is self-limiting (near-identity where u_model is
already close to the trusted/target velocity) and pushes C in the
correct direction (more friction where Albany is running too fast,
less where it is running too slow).

Cells where the update is undefined (e.g. not grounded, mu/N/N_source
invalid, or no valid/positive u_model) are left at C_old unchanged.
"""

import argparse
import json
import os
import shutil

import numpy as np
import xarray as xr

from friction_law_conversion import (
    ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT,
    RC_POWER_EXPONENT,
    SECONDS_PER_YEAR,
    albany_temperature_based_flow_rate,
    compute_named_effective_pressure,
    downs_johnson_effective_pressure,
    effective_pressure4,
    load_albany_ascii_geometry,
    ocean_connection_effective_pressure,
)

# Physics options this script needs to match friction_law_conversion.py's
# original conversion exactly. Values here are friction_law_conversion.py's
# own hardcoded defaults, used only as a last resort when a given option is
# supplied neither on the command line nor via a conversion-config.json
# sidecar file (see --config/parse_args()).
HARD_DEFAULTS = {
    "weertman_q": 0.2,
    "rc_power_exponent": RC_POWER_EXPONENT,
    "glen_n": 3.0,
    "flow_rate_type": "temperature",
    "flow_rate": None,
    "temperature_field": "temperature",
    "effective_pressure_type": "downs-johnson",
    "effective_pressure_input_field": "effectivePressure",
    "min_fraction_overburden": None,
    "pressure_length_scale": None,
    "transition_h_ocean": 25.0,
    "source_effective_pressure_type": "constant",
    "source_effective_pressure": 1.0,
    "source_effective_pressure_field": "effectivePressure",
    "source_min_fraction_overburden": None,
    "source_pressure_length_scale": None,
    "source_transition_h_ocean": 25.0,
    "rho_ice": 910.0,
    "rho_water": 1028.0,
    "gravity": 9.80616,
    "mu_field": "muFriction",
    "lambda_field": "bedRoughnessRC",
    "thickness_field": "thickness",
    "bed_field": "bedTopography",
    "ascii_mesh_dir": ".",
    "velocity_x_field": "uReconstructX",
    "velocity_y_field": "uReconstructY",
    "time_index": 0,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Refine a Regularized Coulomb C (muFriction) field via a "
            "fixed-point (Picard) update against Albany's actual "
            "simulated velocity."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "original_input",
        help=(
            "The *original* MALI initial-condition file -- the same "
            "one originally passed to friction_law_conversion.py -- "
            "used for mu (the source Weertman/Budd friction field), "
            "the source law's N_source, and the original target "
            "velocity/Tau_b_target. Its thickness/bedTopography are "
            "overridden by --ascii-mesh-dir exactly as in "
            "friction_law_conversion.py."
        )
    )
    parser.add_argument(
        "previous_conversion_output",
        help=(
            "friction_law_conversion.py's output file from the pass "
            "whose Lambda/C were actually used for the Albany run "
            "being refined against (i.e. the file supplying C_old "
            "and Lambda; see --mu-field/--lambda-field)."
        )
    )
    parser.add_argument(
        "model_output",
        help=(
            "A MALI-compatible output file containing the *actual* "
            "basal velocity (see --velocity-x-field/"
            "--velocity-y-field) that Albany produced when run with "
            "previous_conversion_output's Lambda/C fields."
        )
    )
    parser.add_argument(
        "output",
        help="Refined output file to write (C updated, Lambda unchanged)."
    )
    parser.add_argument(
        "--config",
        default=None,
        help=(
            "Path to a physics-options JSON config file written by "
            "friction_law_conversion.py (<its output>.config.json), "
            "supplying default values for every physics option below "
            "so they don't need to be re-specified by hand and kept "
            "in sync with the original conversion run. Defaults to "
            "'<previous_conversion_output>.config.json'. Explicit "
            "CLI options always take precedence over the config "
            "file; options given by neither fall back to "
            "friction_law_conversion.py's own hardcoded defaults "
            "(see HARD_DEFAULTS). A warning is printed for any option "
            "where the CLI value overrides a *different* config-file "
            "value, since that usually indicates an unintentional "
            "mismatch."
        )
    )

    parser.add_argument(
        "--relaxation",
        type=float,
        default=1.0,
        help=(
            "Picard step damping: C_refined = C_old * (C_new / "
            "C_old) ** relaxation. 1.0 (default) takes the full "
            "closed-form step; values < 1 damp the update, useful if "
            "repeated iterations oscillate rather than converge."
        )
    )
    parser.add_argument(
        "--bound-factor",
        type=float,
        default=None,
        help=(
            "Optional safety-net bound: any C_refined outside "
            "[C_old / bound_factor, C_old * bound_factor], or non-"
            "finite/non-positive, is reset to C_old instead. "
            "Disabled (no clipping) by default."
        )
    )

    parser.add_argument(
        "--weertman-q", "--q",
        dest="weertman_q",
        type=float,
        default=None,
        help=(
            "Input Weertman/Power-Law sliding exponent qW, used to "
            "compute Tau_b_target (default: friction_law_conversion."
            "py's own default, 0.2, or the config-file value; see "
            "--config). Must match the value used for the original "
            "friction_law_conversion.py run."
        )
    )
    parser.add_argument(
        "--rc-power-exponent",
        type=float,
        default=None,
        help=(
            "Regularized Coulomb law's Power Exponent qR, used to "
            "compute k(u_model) (default: friction_law_conversion."
            "py's own default, 1/3, or the config-file value; see "
            "--config). Must match the value used for the original "
            "friction_law_conversion.py run (see its "
            "--rc-power-exponent)."
        )
    )
    parser.add_argument(
        "--glen-n",
        type=float,
        default=None,
        help=(
            "Glen-law exponent n, used to compute u_c(u_model) "
            "(default: friction_law_conversion.py's own default, 3, "
            "or the config-file value; see --config). Must match the "
            "original friction_law_conversion.py run."
        )
    )
    parser.add_argument(
        "--flow-rate-type",
        choices=["temperature", "constant"],
        default=None,
        help=(
            "How to obtain the Glen flow rate A (default: "
            "friction_law_conversion.py's own default, temperature, "
            "or the config-file value; see --config); must match the "
            "original friction_law_conversion.py run. See that "
            "script's --flow-rate-type for details."
        )
    )
    parser.add_argument(
        "--temperature-field",
        default=None,
        help=(
            "Ice temperature field [K] (default: temperature, or the "
            "config-file value; see --config)"
        )
    )
    parser.add_argument(
        "--flow-rate",
        type=float,
        default=None,
        help=(
            "Constant Glen flow rate A [Pa^-3 s^-1], required when "
            "--flow-rate-type=constant (default: the config-file "
            "value, if any; see --config)."
        )
    )

    parser.add_argument(
        "--effective-pressure-type",
        choices=[
            "downs-johnson", "ocean-connection", "transition", "field",
        ],
        default=None,
        help=(
            "How to compute N (default: friction_law_conversion.py's "
            "own default, downs-johnson, or the config-file value; "
            "see --config); must match the original "
            "friction_law_conversion.py run. See that script's "
            "--effective-pressure-type for details."
        )
    )
    parser.add_argument(
        "--effective-pressure-input-field",
        default=None,
        help=(
            "Field holding a precomputed N [Pa], used only when "
            "--effective-pressure-type=field (default: "
            "effectivePressure, or the config-file value; see "
            "--config). Read from --effective-pressure-file (default: "
            "original_input)."
        )
    )
    parser.add_argument(
        "--effective-pressure-file",
        default=None,
        help=(
            "File to read --effective-pressure-input-field from "
            "when --effective-pressure-type=field (default: "
            "original_input)."
        )
    )
    parser.add_argument(
        "--min-fraction-overburden", type=float, default=None,
        help=(
            "Required for --effective-pressure-type/"
            "--source-effective-pressure-type of downs-johnson or "
            "transition (default: the config-file value, if any; "
            "see --config); see friction_law_conversion.py."
        )
    )
    parser.add_argument(
        "--pressure-length-scale", type=float, default=None,
        help=(
            "Required for --effective-pressure-type/"
            "--source-effective-pressure-type of downs-johnson or "
            "transition (default: the config-file value, if any; "
            "see --config); see friction_law_conversion.py."
        )
    )
    parser.add_argument(
        "--transition-h-ocean", type=float, default=None,
        help=(
            "Only used when --effective-pressure-type=transition "
            "(default: 25.0, or the config-file value; see --config)."
        )
    )

    parser.add_argument(
        "--source-effective-pressure-type",
        choices=[
            "constant", "downs-johnson", "ocean-connection",
            "transition", "field",
        ],
        default=None,
        help=(
            "How to compute N_source for Tau_b_target (default: "
            "friction_law_conversion.py's own default, constant, or "
            "the config-file value; see --config); must match the "
            "original friction_law_conversion.py run."
        )
    )
    parser.add_argument(
        "--source-effective-pressure", type=float, default=None,
        help=(
            "Used when --source-effective-pressure-type=constant "
            "(default: 1.0, or the config-file value; see --config)."
        )
    )
    parser.add_argument(
        "--source-effective-pressure-field",
        default=None,
        help=(
            "Used when --source-effective-pressure-type=field "
            "(default: effectivePressure, or the config-file value; "
            "see --config)."
        )
    )
    parser.add_argument(
        "--source-min-fraction-overburden", type=float, default=None,
        help="Default: the config-file value, if any; see --config."
    )
    parser.add_argument(
        "--source-pressure-length-scale", type=float, default=None,
        help="Default: the config-file value, if any; see --config."
    )
    parser.add_argument(
        "--source-transition-h-ocean", type=float, default=None,
        help=(
            "Default: 25.0, or the config-file value; see --config."
        )
    )

    parser.add_argument("--rho-ice", type=float, default=None)
    parser.add_argument("--rho-water", type=float, default=None)
    parser.add_argument("--gravity", type=float, default=None)

    parser.add_argument("--mu-field", default=None)
    parser.add_argument("--lambda-field", default=None)
    parser.add_argument("--thickness-field", default=None)
    parser.add_argument("--bed-field", default=None)
    parser.add_argument("--ascii-mesh-dir", default=None)
    parser.add_argument(
        "--velocity-x-field", default=None,
        help=(
            "Basal x-velocity field [m s^-1] (default: uReconstructX, "
            "or the config-file value; see --config), read from both "
            "original_input (for the original target speed) and "
            "model_output (for u_model)."
        )
    )
    parser.add_argument(
        "--velocity-y-field", default=None,
        help=(
            "Basal y-velocity field [m s^-1] (default: uReconstructY, "
            "or the config-file value; see --config)"
        )
    )
    parser.add_argument("--time-index", type=int, default=None)

    args = parser.parse_args()

    # -------------------------------------------------------------
    # Merge in the physics-options config file written by
    # friction_law_conversion.py (see HARD_DEFAULTS above and
    # --config's help text): explicit CLI values always win; any
    # option left at None falls back to the config file, then to
    # friction_law_conversion.py's own hardcoded default.
    # -------------------------------------------------------------
    config_path = args.config
    if config_path is None:
        config_path = args.previous_conversion_output + ".config.json"
    config = {}
    if os.path.isfile(config_path):
        with open(config_path) as f:
            config = json.load(f)
        print(f"Loaded physics-options config: {config_path}")
    elif args.config is not None:
        parser.error(f"--config file not found: {config_path!r}")
    else:
        print(
            f"NOTE: no physics-options config file found at "
            f"{config_path!r} (expected alongside "
            "previous_conversion_output, written automatically by a "
            "recent friction_law_conversion.py; see --config). "
            "Falling back to explicit CLI options and/or "
            "friction_law_conversion.py's own hardcoded defaults for "
            "any physics option not given on the command line -- "
            "double check these match the original conversion run."
        )

    for key, hard_default in HARD_DEFAULTS.items():
        cli_value = getattr(args, key)
        config_value = config.get(key)
        if cli_value is not None:
            if config_value is not None and cli_value != config_value:
                print(
                    f"NOTE: --{key.replace('_', '-')}={cli_value!r} "
                    f"on the command line overrides the config-file "
                    f"value {config_value!r} -- make sure this is "
                    "intentional."
                )
            continue
        setattr(args, key, config_value if config_value is not None
                else hard_default)

    if args.flow_rate_type == "constant" and args.flow_rate is None:
        parser.error(
            "--flow-rate is required when --flow-rate-type=constant"
        )
    if args.effective_pressure_type in ("downs-johnson", "transition") and (
        args.min_fraction_overburden is None
        or args.pressure_length_scale is None
    ):
        parser.error(
            "--min-fraction-overburden and --pressure-length-scale "
            "are required when --effective-pressure-type=downs-"
            "johnson or transition"
        )
    if args.source_effective_pressure_type in (
        "downs-johnson", "transition"
    ) and (
        args.source_min_fraction_overburden is None
        or args.source_pressure_length_scale is None
    ):
        parser.error(
            "--source-min-fraction-overburden and --source-pressure-"
            "length-scale are required when "
            "--source-effective-pressure-type=downs-johnson or "
            "transition"
        )

    return args


def cell_field(ds, name, time_index):
    """Extract an nCells field, dropping Time if present."""
    da = ds[name]
    if "Time" in da.dims:
        da = da.isel(Time=time_index)
    values = np.asarray(da.values).squeeze()
    if values.ndim != 1:
        raise ValueError(
            f"{name} must reduce to a 1-D nCells field; got shape "
            f"{values.shape}"
        )
    return values.astype(np.float64)


def basal_cell_field(ds, name, time_index):
    """
    Extract an nCells field from a (Time, nCells, nVertLevels) or
    similar field, taking the last vertical level/interface as an
    approximation of the basal-most value.
    """
    da = ds[name]
    if "Time" in da.dims:
        da = da.isel(Time=time_index)
    vert_dims = [d for d in da.dims if d.lower().startswith("nvert")]
    if vert_dims:
        da = da.isel({vert_dims[0]: -1})
    values = np.asarray(da.values).squeeze()
    if values.ndim != 1:
        raise ValueError(
            f"{name} must reduce to a 1-D nCells field; got shape "
            f"{values.shape}"
        )
    return values.astype(np.float64)


def compute_N(args, H, bed, ds_for_field):
    """Compute N [Pa] the same way friction_law_conversion.py does."""
    if args.effective_pressure_type == "field":
        return cell_field(
            ds_for_field, args.effective_pressure_input_field,
            args.time_index,
        )
    elif args.effective_pressure_type == "ocean-connection":
        return ocean_connection_effective_pressure(
            thickness=H, bed=bed, rho_i=args.rho_ice,
            rho_w=args.rho_water, gravity=args.gravity,
        )
    elif args.effective_pressure_type == "downs-johnson":
        return downs_johnson_effective_pressure(
            thickness=H, bed=bed,
            min_fraction_overburden=args.min_fraction_overburden,
            length_scale=args.pressure_length_scale,
            rho_i=args.rho_ice, rho_w=args.rho_water,
            gravity=args.gravity,
        )
    else:
        return effective_pressure4(
            thickness=H, bed=bed,
            min_fraction_overburden=args.min_fraction_overburden,
            length_scale=args.pressure_length_scale,
            rho_i=args.rho_ice, rho_w=args.rho_water,
            gravity=args.gravity, h_ocean=args.transition_h_ocean,
        )


def main():
    args = parse_args()

    ds_in = xr.open_dataset(args.original_input)
    ds_prev = xr.open_dataset(args.previous_conversion_output)
    ds_model = xr.open_dataset(args.model_output)

    mu = cell_field(ds_in, args.mu_field, args.time_index)
    H = cell_field(ds_in, args.thickness_field, args.time_index)
    bed = cell_field(ds_in, args.bed_field, args.time_index)

    ascii_cell_index, ascii_H, ascii_bed = load_albany_ascii_geometry(
        args.ascii_mesh_dir, n_cells=H.shape[0]
    )
    H[ascii_cell_index] = ascii_H
    bed[ascii_cell_index] = ascii_bed

    grounded = (H > 0.0) & (args.rho_ice * H + args.rho_water * bed > 0.0)

    field_ds = (
        xr.open_dataset(args.effective_pressure_file)
        if args.effective_pressure_file else ds_in
    )
    N = compute_N(args, H, bed, field_ds)
    N_albany = N / ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT

    if args.source_effective_pressure_type == "constant":
        N_source_kpa = np.full_like(H, args.source_effective_pressure)
    elif args.source_effective_pressure_type == "field":
        N_source_kpa = (
            cell_field(
                ds_in, args.source_effective_pressure_field,
                args.time_index,
            ) / ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT
        )
    else:
        N_source_kpa = compute_named_effective_pressure(
            args.source_effective_pressure_type,
            thickness=H, bed=bed, rho_i=args.rho_ice,
            rho_w=args.rho_water, gravity=args.gravity,
            min_fraction_overburden=args.source_min_fraction_overburden,
            pressure_length_scale=args.source_pressure_length_scale,
            transition_h_ocean=args.source_transition_h_ocean,
        ) / ALBANY_EFFECTIVE_PRESSURE_PA_PER_UNIT

    # Original (trusted) target velocity -- used ONLY for the
    # convergence/mismatch diagnostic below, NOT for Tau_b_target
    # itself. See the module docstring: Tau_b_target must be
    # evaluated at Albany's *actual* simulated velocity u_model, not
    # at this original velocity -- otherwise the update is
    # structurally backwards. k(u) is monotonically increasing in u
    # (it saturates at N_albany as u -> infinity), so pinning
    # Tau_b_target at u_target while evaluating k(.) at u_model
    # always pushes C in the wrong direction whenever u_model !=
    # u_target: e.g. if u_model > u_target (too fast), k(u_model) >
    # k(u_target), so C_new = Tau_b_target(u_target) / (C_old *
    # k(u_model)) comes out *smaller* than C_old -- reducing
    # friction and making the too-fast region even faster, which is
    # exactly backwards. Evaluating Tau_b_target at u_model instead
    # (the standard fixed-point/Picard approach: recalibrate C so
    # the RC law agrees with the trusted source law *at whatever
    # velocity Albany is actually producing right now*) fixes this,
    # and self-limits to no-op in already-good regions where u_model
    # ~= u_target.
    uX_target = basal_cell_field(ds_in, args.velocity_x_field, args.time_index)
    uY_target = basal_cell_field(ds_in, args.velocity_y_field, args.time_index)
    u_target = np.sqrt(uX_target ** 2 + uY_target ** 2) * SECONDS_PER_YEAR

    # Actual Albany-simulated velocity from the model run driven by
    # previous_conversion_output's Lambda/C.
    uX_model = basal_cell_field(ds_model, args.velocity_x_field, args.time_index)
    uY_model = basal_cell_field(ds_model, args.velocity_y_field, args.time_index)
    u_model = np.sqrt(uX_model ** 2 + uY_model ** 2) * SECONDS_PER_YEAR

    tau_defined = (
        grounded
        & np.isfinite(N) & (N > 0.0)
        & np.isfinite(mu) & (mu > 0.0)
        & np.isfinite(N_source_kpa) & (N_source_kpa > 0.0)
        & np.isfinite(u_model) & (u_model > 0.0)
    )
    tau_b_target = np.full_like(N, np.nan)
    tau_b_target[tau_defined] = (
        mu[tau_defined] * N_source_kpa[tau_defined]
        * u_model[tau_defined] ** args.weertman_q
    )

    C_old = cell_field(ds_prev, args.mu_field, args.time_index)
    Lambda = cell_field(ds_prev, args.lambda_field, args.time_index)

    # Glen flow rate A, matching friction_law_conversion.py.
    if args.flow_rate_type == "temperature":
        basal_temperature = basal_cell_field(
            ds_in, args.temperature_field, args.time_index
        )
        A = albany_temperature_based_flow_rate(basal_temperature)
    else:
        A = np.full_like(H, args.flow_rate)

    with np.errstate(divide="ignore", invalid="ignore"):
        u_c_model = (
            SECONDS_PER_YEAR * Lambda * A * np.maximum(N, 0.0) ** args.glen_n
        )
        k_model = (
            N_albany * u_model ** args.rc_power_exponent
            / (u_model + u_c_model) ** args.rc_power_exponent
        )

    update_defined = (
        tau_defined
        & np.isfinite(k_model) & (k_model > 0.0)
        & np.isfinite(C_old) & (C_old > 0.0)
    )

    # -------------------------------------------------------------
    # Convergence/mismatch diagnostics (computed early so the
    # self-consistency check below can restrict itself to already-
    # "good" cells).
    # -------------------------------------------------------------
    vel_mismatch = np.full_like(H, np.nan)
    mismatch_defined = (
        grounded
        & np.isfinite(u_model) & (u_model > 0.0)
        & np.isfinite(u_target) & (u_target > 0.0)
    )
    vel_mismatch[mismatch_defined] = (
        (u_model[mismatch_defined] - u_target[mismatch_defined])
        / u_target[mismatch_defined]
    )

    # -------------------------------------------------------------
    # Self-consistency check (parameter-mismatch diagnostic): using
    # THIS SCRIPT's own N/N_source/u_c and the ORIGINAL target
    # velocity u_target (not u_model), recompute what C_old *should*
    # be, per the very same formula friction_law_conversion.py's
    # solve_transition_velocity() used to derive it in the first
    # place. If all physics options passed to this script (--weertman
    # -q, --glen-n, --rc-power-exponent, --effective-pressure-type
    # and related thresholds, --source-effective-pressure-type and
    # related thresholds, --flow-rate-type/--flow-rate, rho/gravity
    # values, --ascii-mesh-dir) match those used for the ORIGINAL
    # friction_law_conversion.py run exactly, C_check should equal
    # C_old almost exactly (to floating-point precision) everywhere
    # C_old came from that same closed-form solve. A large mismatch
    # here -- especially in cells where u_model ~= u_target (i.e.
    # already "good" regions) -- means at least one physics option
    # differs between the two runs, NOT that the Picard update itself
    # is misbehaving; a mismatch confined to cells that are *not*
    # already "good" is expected/uninformative, since C_old there may
    # have come from a different code path (e.g.
    # fit_coulomb_C_fast_region()'s scalar fit for
    # --method=stress-match-fit) that this check does not reproduce.
    # -------------------------------------------------------------
    k_target = np.full_like(N, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        k_target[tau_defined] = (
            N_albany[tau_defined]
            * u_target[tau_defined] ** args.rc_power_exponent
            / (u_target[tau_defined] + u_c_model[tau_defined])
            ** args.rc_power_exponent
        )
    tau_b_at_target = np.full_like(N, np.nan)
    tau_b_at_target[tau_defined] = (
        mu[tau_defined] * N_source_kpa[tau_defined]
        * u_target[tau_defined] ** args.weertman_q
    )
    C_check = np.full_like(N, np.nan)
    check_defined = (
        tau_defined & np.isfinite(k_target) & (k_target > 0.0)
        & np.isfinite(C_old) & (C_old > 0.0)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        C_check[check_defined] = (
            tau_b_at_target[check_defined] / k_target[check_defined]
        )
    consistency_residual = np.full_like(N, np.nan)
    consistency_residual[check_defined] = (
        np.abs(C_check[check_defined] - C_old[check_defined])
        / C_old[check_defined]
    )
    already_good = (
        mismatch_defined & check_defined & (np.abs(vel_mismatch) < 0.05)
    )
    if np.count_nonzero(already_good) > 0:
        print(
            "Self-consistency check (C_check vs. C_old, using "
            "u_target in place of u_model -- should match almost "
            "exactly if all physics options match the original "
            "friction_law_conversion.py run) for cells where "
            "|model - target| / target < 5% : "
            f"mean={np.nanmean(consistency_residual[already_good]):.6e}, "
            f"95th pct={np.nanpercentile(consistency_residual[already_good], 95):.6e}, "
            f"max={np.nanmax(consistency_residual[already_good]):.6e} "
            "(large values here indicate a physics-option mismatch "
            "between this call and the original conversion run, not "
            "a problem with the Picard update itself)"
        )

    # Direct closed-form C that makes the RC law reproduce
    # Tau_b_target *at u_model* -- i.e. exactly
    # solve_transition_velocity()'s C formula, but re-solved using
    # Albany's actual simulated velocity in place of the original
    # conversion's speed. This is a full (unrelaxed) Picard step; see
    # --relaxation for damping it.
    C_new_raw = np.full_like(C_old, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        C_new_raw[update_defined] = (
            tau_b_target[update_defined] / k_model[update_defined]
        )

    C_new = C_old.copy()
    ratio = np.full_like(C_old, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio[update_defined] = (
            C_new_raw[update_defined] / C_old[update_defined]
        )
    C_new[update_defined] = (
        C_old[update_defined] * ratio[update_defined] ** args.relaxation
    )

    n_reset = 0
    if args.bound_factor is not None:
        lower = C_old / args.bound_factor
        upper = C_old * args.bound_factor
        reset = (
            update_defined
            & (
                ~np.isfinite(C_new)
                | (C_new <= 0.0)
                | (C_new < lower)
                | (C_new > upper)
            )
        )
        C_new[reset] = C_old[reset]
        n_reset = int(np.count_nonzero(reset))

    n_grounded = int(np.count_nonzero(grounded))
    n_updated = int(np.count_nonzero(update_defined))
    print(
        f"Relative velocity mismatch ((model - target) / target), "
        f"grounded : mean(|.|)={np.nanmean(np.abs(vel_mismatch[grounded])):.6e}, "
        f"95th pct={np.nanpercentile(np.abs(vel_mismatch[grounded]), 95):.6e}, "
        f"max={np.nanmax(np.abs(vel_mismatch[grounded])):.6e}"
    )
    print(
        f"C refined for {n_updated} / {n_grounded} grounded cells "
        f"({n_reset} reset to C_old by --bound-factor)"
    )
    print(
        "C range before/after refinement (grounded, updated cells) : "
        f"{np.nanmin(C_old[update_defined]):.6e} -- "
        f"{np.nanmax(C_old[update_defined]):.6e} (before) -> "
        f"{np.nanmin(C_new[update_defined]):.6e} -- "
        f"{np.nanmax(C_new[update_defined]):.6e} (after)"
        if n_updated > 0 else
        "C range before/after refinement : n/a (no cells updated)"
    )
    if n_updated > 0:
        print(
            "Picard ratio (C_new/C_old before relaxation), grounded "
            "updated cells : "
            f"min={np.nanmin(ratio[update_defined]):.6e}, "
            f"median={np.nanmedian(ratio[update_defined]):.6e}, "
            f"max={np.nanmax(ratio[update_defined]):.6e}"
        )

    # -------------------------------------------------------------
    # Write output: a copy of the *original* input, with C
    # (--mu-field) refined and Lambda (--lambda-field) carried over
    # unchanged from previous_conversion_output -- matching
    # friction_law_conversion.py's own output convention.
    # -------------------------------------------------------------
    ds_in.close()
    ds_prev.close()
    ds_model.close()
    shutil.copy2(args.original_input, args.output)

    out = xr.open_dataset(args.output).load()
    ncell_dim = None
    for dim in out[args.mu_field].dims:
        if dim.lower() == "ncells":
            ncell_dim = dim
            break
    if ncell_dim is None:
        ncell_dim = "nCells"

    out[args.mu_field] = xr.DataArray(
        C_new, dims=(ncell_dim,),
        attrs={
            "long_name": (
                "Albany regularized-Coulomb coefficient C (Mu), "
                "refined via a fixed-point Picard update against "
                "Albany's own actual simulated velocity -- see "
                "refine_regularized_coulomb_mu.py"
            ),
            "units": "1",
        },
    )
    out[args.lambda_field] = xr.DataArray(
        Lambda, dims=(ncell_dim,),
        attrs={
            "long_name": (
                "Albany regularized-Coulomb bed roughness Lambda, "
                "carried over unchanged from "
                f"{args.previous_conversion_output!r}"
            ),
            "units": "m",
        },
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        refinement_ratio = np.where(update_defined, C_new / C_old, np.nan)
    out["muFrictionRefinementRatio"] = xr.DataArray(
        refinement_ratio,
        dims=(ncell_dim,),
        attrs={
            "long_name": (
                "C_refined / C_old from this fixed-point refinement "
                "pass (NaN where not updated)"
            ),
            "units": "1",
        },
    )
    out["velocityMismatchRelative"] = xr.DataArray(
        vel_mismatch, dims=(ncell_dim,),
        attrs={
            "long_name": (
                "(model_output velocity - original_input target "
                "velocity) / target velocity, at the model_output "
                "velocity that resulted from previous_conversion_"
                "output's Lambda/C"
            ),
            "units": "1",
        },
    )
    # xarray cannot safely overwrite an open source file, so use a temp
    # file, matching friction_law_conversion.py's own convention.
    tmp = args.output + ".tmp"
    out.to_netcdf(tmp)
    out.close()
    shutil.move(tmp, args.output)

    print(f"Wrote refined IC: {args.output}")


if __name__ == "__main__":
    main()
