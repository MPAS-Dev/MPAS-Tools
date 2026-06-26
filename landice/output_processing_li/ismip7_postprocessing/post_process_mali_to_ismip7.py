#!/usr/bin/env python

"""
This script processes MALI simulation outputs into the required format
for submission to ISMIP7.
There are 3 flavors of variables that can be processed:
1D, 2D state, and 2D flux variables.
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

Example usage from testing on initial ISMIP7 submission data:
python post_process_mali_to_ismip7.py \
    -o /pscratch/sd/h/hoffman2/ISMIP7-postprocessing-June30-deadline/AIS/test \
    -g "/pscratch/sd/t/trhille/ISMIP7/AIS_runs/no_slm_20260616/landice/ismip7_run/ismip7_ais/historical_CESM2-WACCM/output/globalStats_*.nc" \
    --icesheet AIS \
    -e C001 \
    --res 04 \
    -m /pscratch/sd/t/trhille/ISMIP7/AIS_runs/no_slm_20260616/landice/ismip7_run/ismip7_ais/historical_CESM2-WACCM/relaxed_10yrs_4km.nc \
    -s "/pscratch/sd/t/trhille/ISMIP7/AIS_runs/no_slm_20260616/landice/ismip7_run/ismip7_ais/historical_CESM2-WACCM/output/output_state_*.nc" \
    --ismip7_grid_file /pscratch/sd/h/hoffman2/ISMIP7-postprocessing-June30-deadline/AIS/misc/af2_AIS_04000m_v1.nc \
    -f "/pscratch/sd/t/trhille/ISMIP7/AIS_runs/no_slm_20260616/landice/ismip7_run/ismip7_ais/historical_CESM2-WACCM/output/output_state_*.nc"

And here is an example of only processing 2d state and flux after the mapping file has been created:
python post_process_mali_to_ismip7.py \
    -o /pscratch/sd/h/hoffman2/ISMIP7-postprocessing-June30-deadline/AIS/test \
    --icesheet AIS \
    -e C001 \
    --res 04 \
    -m /pscratch/sd/t/trhille/ISMIP7/AIS_runs/no_slm_20260616/landice/ismip7_run/ismip7_ais/historical_CESM2-WACCM/relaxed_10yrs_4km.nc \
    -s "/pscratch/sd/t/trhille/ISMIP7/AIS_runs/no_slm_20260616/landice/ismip7_run/ismip7_ais/historical_CESM2-WACCM/output/output_state_*.nc" \
    -f "/pscratch/sd/t/trhille/ISMIP7/AIS_runs/no_slm_20260616/landice/ismip7_run/ismip7_ais/historical_CESM2-WACCM/output/output_state_*.nc" \
    --ismip7_grid_file /pscratch/sd/h/hoffman2/ISMIP7-postprocessing-June30-deadline/AIS/misc/af2_AIS_04000m_v1.nc \
    --reuse_mapping_file mapping_mali_to_ismip7.conserve.20260626T130744.nc

"""

import argparse
import glob
import os
from datetime import datetime

from grid_and_mapping import (
    build_mapping_file,
    check_exp_name,
    check_ismip7_grid_file,
    check_res,
)
from process_1d_variables_ismip7 import (
    check_global_stats_files,
    generate_output_1d_vars,
)
from process_flux_variables_ismip7 import process_flux_pipeline
from process_state_variables_ismip7 import process_state_pipeline

DEFAULT_AUTHORS = 'Matthew Hoffman, Trevor Hillebrand, Holly Kyeore Han'
DEFAULT_GROUP = 'Los Alamos National Laboratory, Department of Energy'
DEFAULT_MODEL = 'MALI (MPAS-Albany Land Ice model)'
DEFAULT_GROUP_NICKNAME = 'DOE'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-e",
        "--exp_name",
        dest="exp",
        required=True,
        help="ISMIP7 experiment name (e.g., exp05)",
    )
    parser.add_argument(
        "-s",
        "--input_state_pattern",
        dest="input_state_pattern",
        required=False,
        help=(
            "glob pattern matching one or more MALI state output files. "
            "Note: wildcard paths should be quoted to avoid shell expansion."
        ),
    )
    parser.add_argument(
        "-f",
        "--input_flux_pattern",
        dest="input_flux_pattern",
        required=False,
        help=(
            "glob pattern matching one or more MALI flux output files. "
            "Note: wildcard paths should be quoted to avoid shell expansion."
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
        "-g",
        "--global_stats_pattern",
        dest="global_stats_pattern",
        required=False,
        help=(
            "glob pattern matching one or more globalStats.nc files "
            "(e.g. 'globalStats_*.nc'). "
            "Note: wildcard paths should be quoted to avoid shell expansion."
        ),
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
        "--ismip7_grid_file",
        dest="ismip7_grid_file",
        required=True,
        help="Input ismip7 mesh file.",
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
        help="resolution of the ismip7 grid, in kilometers: 16, 8, 4, 2, 1",
    )
    parser.add_argument(
        "--icesheet",
        dest="icesheet",
        required=True,
        choices=['AIS', 'GIS'],
        help="ice sheet domain: 'AIS' (Antarctica) or 'GIS' (Greenland)",
    )
    parser.add_argument(
        "--authors",
        dest="authors",
        required=False,
        default=DEFAULT_AUTHORS,
        help=(
            "author string for output file metadata "
            f"(default: '{DEFAULT_AUTHORS}')"
        ),
    )
    parser.add_argument(
        "--group",
        dest="group",
        required=False,
        default=DEFAULT_GROUP,
        help=(
            "group/institution string for output file metadata "
            f"(default: '{DEFAULT_GROUP}')"
        ),
    )
    parser.add_argument(
        "--group_nickname",
        dest="group_nickname",
        required=False,
        default=DEFAULT_GROUP_NICKNAME,
        help=(
            "short group nickname used in output filenames "
            f"(default: '{DEFAULT_GROUP_NICKNAME}')"
        ),
    )
    args = parser.parse_args()

    check_exp_name(args.exp)
    check_res(args.res_ismip7_grid)

    metadata = {
        'exp': args.exp,
        'icesheet': args.icesheet,
        'authors': args.authors,
        'group': args.group,
        'group_nickname': args.group_nickname,
        'model': DEFAULT_MODEL,
        'date': datetime.now().strftime("%d-%b-%Y"),
    }

    print("\n---Processing remapping file---")
    # Only do remapping steps if we have 2d files to process
    if (
            args.input_state_pattern is not None or
            args.input_flux_pattern is not None):

        check_ismip7_grid_file(args.ismip7_grid_file, args.res_ismip7_grid)

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
            mapping_file = (
                f"mapping_mali_to_ismip7.{method_remap}.{created_at}.nc"
            )

            print("Creating new mapping file. "
                  f"Mapping method used: {method_remap}")
            build_mapping_file(
                args.input_file_mesh,
                mapping_file,
                args.res_ismip7_grid,
                args.ismip7_grid_file,
                method_remap,
            )

    print("---Processing remapping file complete---\n")

    output_path = args.output_path
    print(f"Using output path: {output_path}")
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    if args.global_stats_pattern is None:
        print("--- No global stats pattern provided; "
              "skipping 1D variable processing.")
    else:
        print("\n---Processing global stats file(s)---")
        global_stats_files = sorted(glob.glob(args.global_stats_pattern))
        check_global_stats_files(global_stats_files)
        generate_output_1d_vars(global_stats_files, output_path, metadata)
        print("---Processing global stats file(s) complete---\n")

    if args.input_state_pattern is None:
        print("--- MALI state pattern is not provided, "
              "thus it will not be processed.")
    else:
        print("\n---Processing state file(s)---")
        state_files = sorted(glob.glob(args.input_state_pattern))
        process_state_pipeline(
            state_files,
            mapping_file,
            args.ismip7_grid_file,
            output_path,
            metadata,
        )
        print("---Processing state file(s) complete---\n")

    if args.input_flux_pattern is None:
        print("--- MALI flux pattern is not provided, "
              "thus it will not be processed.")
    else:
        print("\n---Processing flux file(s)---")
        flux_files = sorted(glob.glob(args.input_flux_pattern))
        process_flux_pipeline(
            flux_files,
            mapping_file,
            args.ismip7_grid_file,
            output_path,
            metadata,
        )
        print("---Processing flux file(s) complete---\n")

    print("---All processing complete---")


if __name__ == "__main__":
    main()
