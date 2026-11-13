#!/usr/bin/env python

"""
This script processes MALI simulation outputs into the required format
for submission to ISMIP7.
There are 4 flavors of variables that can be processed:
1D, 2D state, 2D flux, and restart-derived (litempavg) variables.
Any combination can be specified.
Processing assumes simulations were run with the output streams defined
by the ismip7_run compass test case.
Multifile output file sets should be specified with a glob pattern.
Note: wildcard paths should be quoted to avoid shell expansion.

An ISMIP submimssion grid resolution and grid file need to be specified.
A mapping file will be created unless an existing one is specified.

More info at:
https://github.com/ismip/ISM_SimulationChecker/blob/main/conventions/ISMIP7_variable_request.csv
https://www.ismip.org/research/ismip7

Example usage from testing on GrIS production ensemble:
python post_process_mali_to_ismip7.py -o /pscratch/sd/h/hoffman2/ISMIP7/GIS_PRODUCTION_20261108/landice/ismip7_run/ismip7_gris/postprocessing --icesheet GrIS --res 01 --group DOE -m /global/cfs/cdirs/fanssie/MALI_projects/ISMIP7/GIS/initial_condition/20261029/GIS_1to10km_r03_20260925_after_8yr-iter-relax_basicFields2.nc -s "output_2d_*.nc" -f "output_2d_*.nc" -i /pscratch/sd/h/hoffman2/ISMIP7/GIS_PRODUCTION_20261108/landice/ismip7_run/ismip7_gris/ocx/output

The experiment is determined from the name of the run directory (the parent of
the --mali_output_dir directory), which must be '<experiment>_<ESM>', e.g.
'historical_CESM2-WACCM' or 'ssp126_MRI-ESM2-0', or just 'ocx'.  The name is
matched to experiments.csv, which gives the set counter (e.g. C001) used in the
output file names and output subdirectory.  Use --run_name to override the
directory name.
"""

import argparse
import csv
import glob
import os
from datetime import datetime

from grid_and_mapping import (
    build_mapping_file,
    check_res,
    create_ismip7_grid_file,
    CRS_DICT,
    get_time_range,
)
from process_1d_variables_ismip7 import (
    check_global_stats_files,
    generate_output_1d_vars,
)
from process_flux_variables_ismip7 import process_flux_pipeline
from process_restart_variables_ismip7 import process_restart_pipeline
from process_state_variables_ismip7 import process_state_pipeline

DEFAULT_MODEL = 'MALI (MPAS-Albany Land Ice model)'
DEFAULT_ISM_ID = 'MALI7'

GROUP_METADATA = {
    'DOE': {
        'contact_name': 'Trevor Hillebrand, Matthew Hoffman',
        'contact_email': 'trhille@lanl.gov, mhoffman@lanl.gov',
        'group': 'Los Alamos National Laboratory, U.S. Department of Energy',
        'group_nickname': 'DOE',
    },
    'Arete': {
        'contact_name': 'Kyeore Han, Colin Meyer',
        'contact_email': 'hollyhan4@gmail.com, Colin.R.Meyer@dartmouth.edu',
        'group': 'Arete Glacier Initiative',
        'group_nickname': 'ARETE',
    },
}
DEFAULT_ISM_MEMBER_ID = 'm001'
DEFAULT_FORCING_MEMBER_ID = 'f001'
EXPERIMENTS_CSV = os.path.join(os.path.dirname(__file__), 'experiments.csv')


# run directory experiment names that differ from experiments.csv
EXPERIMENT_ALIASES = {'ctrl': 'ctrl2015'}


def _read_experiments():
    with open(EXPERIMENTS_CSV, newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def find_experiment(run_name):
    """
    Match a run name of the form '<experiment>_<ESM>' (or just
    '<experiment>' for experiments without an ESM, e.g. 'ocx') to a row of
    experiments.csv.  Matching is case insensitive.
    """
    rows = _read_experiments()
    valid = [r['experiment_id'] if r['esm_id'].strip() == '-'
             else f"{r['experiment_id']}_{r['esm_id']}" for r in rows]
    if '_' in run_name:
        exp_token, esm_token = run_name.split('_', 1)
    else:
        exp_token, esm_token = run_name, '-'
    exp_token = exp_token.lower()
    exp_token = EXPERIMENT_ALIASES.get(exp_token, exp_token)
    matches = [r for r in rows
               if r['experiment_id'].strip().lower() == exp_token and
               r['esm_id'].strip().lower() == esm_token.lower()]
    if len(matches) != 1:
        raise ValueError(
            f"Could not match run name '{run_name}' to a unique experiment "
            f"in '{EXPERIMENTS_CSV}'. Run directories must be named "
            "'<experiment>_<ESM>' (or just '<experiment>' for ocx), e.g.: "
            f"{', '.join(valid)}. Use --run_name to override the directory "
            "name.")
    return matches[0]


def _run_name_from_output_dir(mali_output_dir):
    """Name of the run directory: the parent of the 'output' directory."""
    path = os.path.abspath(mali_output_dir.rstrip('/'))
    if os.path.basename(path) == 'output':
        path = os.path.dirname(path)
    return os.path.basename(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run_name",
        dest="run_name",
        required=False,
        help=(
            "Run name '<experiment>_<ESM>' (e.g., historical_CESM2-WACCM) "
            "used to look up the experiment in experiments.csv. Default: "
            "name of the run directory containing --mali_output_dir."
        ),
    )
    parser.add_argument(
        "-i",
        "--mali_output_dir",
        dest="mali_output_dir",
        required=True,
        help="Path to MALI output directory containing output NetCDF files.",
    )
    parser.add_argument(
        "-g",
        "--global_stats_pattern",
        dest="global_stats_pattern",
        required=False,
        default='globalStats_*.nc',
        help=(
            "Filename pattern within --mali_output_dir for one or more "
            "globalStats files (default: 'globalStats_*.nc')."
        ),
    )
    parser.add_argument(
        "-s",
        "--input_state_pattern",
        dest="input_state_pattern",
        required=False,
        default='output_2d_*.nc',
        help=(
            "Filename pattern within --mali_output_dir for one or more "
            "MALI 2d state output files (default: 'output_state_*.nc')."
        ),
    )
    parser.add_argument(
        "-f",
        "--input_flux_pattern",
        dest="input_flux_pattern",
        required=False,
        default='output_2d_*.nc',
        help=(
            "Filename pattern within --mali_output_dir for one or more "
            "MALI 2d flux output files (default: 'output_flux_*.nc')."
        ),
    )
    parser.add_argument(
        "--process_sets",
        dest="process_sets",
        nargs='+',
        choices=['1d', 'state', 'flux', 'restart'],
        default=['1d', 'state', 'flux', 'restart'],
        help=(
            "Which variable sets to process. Provide one or more of: "
            "'1d', 'state', 'flux', 'restart'. Default processes all. "
            "'restart' (litempavg) needs restart files; omit it if they "
            "are not available."
        ),
    )
    parser.add_argument(
        "--restart_dir",
        dest="restart_dir",
        required=False,
        help=(
            "Directory containing MALI restart files (default: parent "
            "directory of --mali_output_dir)."
        ),
    )
    parser.add_argument(
        "--input_restart_pattern",
        dest="input_restart_pattern",
        required=False,
        default='rst.*-01-01.nc',
        help=(
            "Filename pattern within --restart_dir for Jan. 1 restart "
            "files (default: 'rst.*-01-01.nc')."
        ),
    )
    parser.add_argument(
        "--fix_ligroundf_thickness_division",
        dest="fix_ligroundf",
        action="store_true",
        default=False,
        help=(
            "Fallback for runs made with an older MALI version whose "
            "avgGroundingLineFlux was divided by the grounded-cell thickness "
            "(huge values at thin ice). Undoes the division, approximately: "
            "flux * thickness * dvEdge / areaCell, using the mean of the "
            "thickness in the flux files and in any restart files within "
            "each output period. Cells beyond +/-10 kg m-2 s-1 (not "
            "recoverable) are set to 0. Requires --input_mesh and applies to "
            "the 'flux' set. Restarts are looked up in --restart_dir with "
            "--ligroundf_restart_pattern (optional). Off by default."
        ),
    )
    parser.add_argument(
        "--ligroundf_restart_pattern",
        dest="ligroundf_restart_pattern",
        required=False,
        default='rst.*.nc',
        help=(
            "Filename pattern within --restart_dir for the restart files "
            "whose thickness is used by --fix_ligroundf_thickness_division "
            "(default: 'rst.*.nc', i.e. all restarts)."
        ),
    )
    parser.add_argument(
        "-m",
        "--input_mesh",
        dest="input_file_mesh",
        required=False,
        help="MALI file with mesh information",
    )
    parser.add_argument(
        "-o",
        "--output_path",
        dest="output_path",
        required=True,
        help="path to which the final output files will be saved",
    )
    parser.add_argument(
        "--reuse_mapping_file",
        dest="reuse_mapping_file",
        required=False,
        help="existing mapping file name to reuse",
    )
    parser.add_argument(
        "--method",
        dest="method_remap",
        default="conserve",
        required=False,
        help="mapping method. Default='conserve'",
    )
    parser.add_argument(
        "--res",
        dest="res_ismip7_grid",
        required=True,
        help=(
            "Resolution of the ISMIP7 output grid in km. "
            "AIS valid values: 2, 4, 8, 16. "
            "GrIS valid values: 1, 2, 4, 8, 16."
        ),
    )
    parser.add_argument(
        "--icesheet",
        dest="icesheet",
        required=True,
        choices=['AIS', 'GrIS'],
        help="ice sheet domain: 'AIS' (Antarctica) or 'GrIS' (Greenland)",
    )
    parser.add_argument(
        "--group",
        dest="group",
        required=True,
        choices=list(GROUP_METADATA.keys()),
        help="Modelling group. Sets authors, institution, and nickname for "
             f"output filenames/metadata. Choices: {list(GROUP_METADATA.keys())}",
    )
    parser.add_argument(
        "--reference_geoid",
        dest="reference_geoid",
        default="EIGEN-6C4",
        help=(
            "Value of the reference_geoid attribute written on deltag: "
            "'EIGEN-6C4' (BedMachine, default), 'GL04C' (Bedmap), or "
            "'custom, defined by ref_geoid variable'."
        ),
    )
    parser.add_argument(
        "--ism_member_id",
        dest="ism_member_id",
        required=False,
        default=DEFAULT_ISM_MEMBER_ID,
        help=(
            "ISM configuration/member identifier for output filenames "
            f"(default: '{DEFAULT_ISM_MEMBER_ID}')"
        ),
    )
    parser.add_argument(
        "--forcing_member_id",
        dest="forcing_member_id",
        required=False,
        default=DEFAULT_FORCING_MEMBER_ID,
        help=(
            "Forcing configuration/member identifier for output filenames "
            f"(default: '{DEFAULT_FORCING_MEMBER_ID}')"
        ),
    )
    args = parser.parse_args()

    process_1d = '1d' in args.process_sets
    process_state = 'state' in args.process_sets
    process_flux = 'flux' in args.process_sets
    process_restart = 'restart' in args.process_sets

    if not os.path.isdir(args.mali_output_dir):
        raise FileNotFoundError(
            f"MALI output directory not found: {args.mali_output_dir}"
        )

    global_stats_glob = os.path.join(
        args.mali_output_dir, args.global_stats_pattern)
    state_glob = os.path.join(args.mali_output_dir, args.input_state_pattern)
    flux_glob = os.path.join(args.mali_output_dir, args.input_flux_pattern)

    global_stats_files = sorted(glob.glob(global_stats_glob))
    check_global_stats_files(global_stats_files)

    check_res(args.res_ismip7_grid, args.icesheet)

    ismip7_grid_file = (
        f"ismip7_grid_{args.icesheet}_{args.res_ismip7_grid}km.nc"
    )
    create_ismip7_grid_file(
        args.icesheet, args.res_ismip7_grid, ismip7_grid_file
    )

    run_name = args.run_name or _run_name_from_output_dir(
        args.mali_output_dir)
    experiment_info = find_experiment(run_name)
    set_counter = experiment_info['set_counter'].strip()
    print(f"Run name '{run_name}' matches experiment {set_counter}: "
          f"experiment_id={experiment_info['experiment_id'].strip()}, "
          f"esm_id={experiment_info['esm_id'].strip()}")
    experiment_id = experiment_info['experiment_id'].strip().lower()
    esm_id = experiment_info['esm_id'].strip()

    group_meta = GROUP_METADATA[args.group]

    metadata = {
        'set_id': 'CORE',
        'set_counter': set_counter,
        'domain_id': args.icesheet,
        'source_id': group_meta['group_nickname'],
        'ism_id': DEFAULT_ISM_ID,
        'ism_member_id': args.ism_member_id,
        'esm_id': esm_id,
        'forcing_member_id': args.forcing_member_id,
        'experiment_id': experiment_id,
        'contact_name': group_meta['contact_name'],
        'contact_email': group_meta['contact_email'],
        'group': group_meta['group'],
        'group_nickname': group_meta['group_nickname'],
        'model': DEFAULT_MODEL,
        'date': datetime.now().strftime("%d-%b-%Y"),
        'crs': CRS_DICT[args.icesheet],
    }

    # Compute one time_range from globalStats and use for all outputs.
    metadata['time_range'] = get_time_range(global_stats_files)
    print(f"Time range: {metadata['time_range']}")

    print("\n---Processing remapping file---")
    # Only do remapping steps if we have 2d files to process
    if process_state or process_flux or process_restart:

        method_remap = args.method_remap
        if args.reuse_mapping_file is not None:
            if not os.path.exists(args.reuse_mapping_file):
                raise FileNotFoundError(
                    "Mapping file to reuse not found: "
                    f"{args.reuse_mapping_file}"
                )
            print(f"Reusing existing mapping file: {args.reuse_mapping_file}")
            mapping_file = args.reuse_mapping_file
        else:
            if args.input_file_mesh is None:
                raise ValueError(
                    "--input_mesh is required when creating "
                    "a new mapping file."
                )

            created_at = datetime.now().strftime("%Y%m%dT%H%M%S")
            mesh_name = os.path.splitext(
                os.path.basename(args.input_file_mesh)
            )[0].replace(' ', '-')
            mapping_file = (
                f"mapping_mali_to_ismip7.{args.icesheet}.{mesh_name}."
                f"{method_remap}.{created_at}.nc"
            )

            print("Creating new mapping file. "
                  f"Mapping method used: {method_remap}")
            print(f"Mapping file name: {mapping_file}")
            build_mapping_file(
                args.input_file_mesh,
                mapping_file,
                args.res_ismip7_grid,
                icesheet=args.icesheet,
                ismip7_grid_file=ismip7_grid_file,
                method_remap=method_remap,
            )

    print("---Processing remapping file complete---\n")

    output_path = os.path.join(args.output_path, set_counter)
    print(f"Using output path: {output_path}")
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    if process_1d:
        print("\n---Processing global stats file(s)---")
        generate_output_1d_vars(global_stats_files, output_path, metadata)
        print("---Processing global stats file(s) complete---\n")
    else:
        print("--- Skipping 1D variable processing (not selected).")

    if process_state:
        print("\n---Processing state file(s)---")
        state_files = sorted(glob.glob(state_glob))
        process_state_pipeline(
            state_files,
            mapping_file,
            ismip7_grid_file,
            output_path,
            metadata,
            reference_geoid=args.reference_geoid,
        )
        print("---Processing state file(s) complete---\n")
    else:
        print("--- Skipping state variable processing (not selected).")

    if process_flux:
        print("\n---Processing flux file(s)---")
        flux_files = sorted(glob.glob(flux_glob))
        gl_restart_files = None
        if args.fix_ligroundf:
            if not args.input_file_mesh:
                raise ValueError(
                    "--input_mesh is required with "
                    "--fix_ligroundf_thickness_division.")
            restart_dir = args.restart_dir or os.path.dirname(
                os.path.abspath(args.mali_output_dir.rstrip('/')))
            gl_restart_files = sorted(glob.glob(
                os.path.join(restart_dir, args.ligroundf_restart_pattern)))
            if not gl_restart_files:
                print("No restart files match "
                      f"{restart_dir}/{args.ligroundf_restart_pattern}; "
                      "using only the thickness in the flux files.")
            print("NOTE: ligroundf is corrected approximately for the old "
                  "division by thickness.")
        process_flux_pipeline(
            flux_files,
            mapping_file,
            ismip7_grid_file,
            output_path,
            metadata,
            input_file_mesh=args.input_file_mesh,
            gl_restart_files=gl_restart_files,
            fix_gl_flux=args.fix_ligroundf,
        )
        print("---Processing flux file(s) complete---\n")
    else:
        print("--- Skipping flux variable processing (not selected).")

    if process_restart:
        print("\n---Processing restart file(s)---")
        if args.input_file_mesh is None:
            raise ValueError(
                "--input_mesh is required for restart processing "
                "(layerThicknessFractions).")
        restart_dir = args.restart_dir or os.path.dirname(
            os.path.abspath(args.mali_output_dir.rstrip('/')))
        restart_files = sorted(glob.glob(
            os.path.join(restart_dir, args.input_restart_pattern)))
        process_restart_pipeline(
            restart_files,
            args.input_file_mesh,
            mapping_file,
            ismip7_grid_file,
            output_path,
            metadata,
        )
        print("---Processing restart file(s) complete---\n")
    else:
        print("--- Skipping restart variable processing (not selected).")

    print("---All processing complete---")


if __name__ == "__main__":
    main()
