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
"""

import argparse
from subprocess import check_call
import os
import glob
from datetime import datetime
from grid_and_mapping import build_mapping_file, check_ismip7_grid_file, \
    check_exp_name, check_res
from process_state_variables_ismip7 import generate_output_2d_state_vars, \
     process_state_vars
from process_1d_variables_ismip7 import generate_output_1d_vars, \
    check_global_stats_files
from process_flux_variables_ismip7 import generate_output_2d_flux_vars

DEFAULT_AUTHORS = 'Matthew Hoffman, Trevor Hillebrand, Holly Kyeore Han'
DEFAULT_GROUP = 'Los Alamos National Laboratory, Department of Energy'
DEFAULT_MODEL = 'MALI (MPAS-Albany Land Ice model)'
DEFAULT_GROUP_NICKNAME = 'DOE'

def main():
    parser = argparse.ArgumentParser(
                        description=__doc__)
    parser.add_argument("-e", "--exp_name", dest="exp",
                        required=True,
                        help="ISMIP7 experiment name (e.g., exp05")
    parser.add_argument("-s", "--input_state", dest="input_file_state",
                        required=False, help="mpas output state variables")
    parser.add_argument("-f", "--input_flux", dest="input_file_flux",
                        required=False, help="mpas output flux variables")
    parser.add_argument("-m", "--input_mesh", dest="input_file_mesh",
                        required=False, help="MALI file with mesh information")
    parser.add_argument("-g", "--global_stats_pattern", dest="global_stats_pattern",
                        required=False,
                        help="glob pattern matching one or more globalStats.nc "
                             "files (e.g. 'globalStats_*.nc')."
                             "Note: wildcard paths should be quoted to avoid shell expansion.")
    parser.add_argument("-o", "--output_path", dest="output_path",
                        required=True,
                        help="path to which the final output files"
                             " will be saved")
    parser.add_argument("--reuse_mapping_file", dest="reuse_mapping_file",
                        required=False,
                        help="existing mapping file name to reuse")
    parser.add_argument("--ismip7_grid_file", dest="ismip7_grid_file",
                        help="Input ismip7 mesh file.")
    parser.add_argument("--method", dest="method_remap", default="conserve",
                        required=False,
                        help="mapping method. Default='conserve'")
    parser.add_argument("--res", dest="res_ismip7_grid",
                        required=True,
                        help="resolution of the ismip7 grid, in kilometers: 16, 8, 4, 2, 1")
    parser.add_argument("--icesheet", dest="icesheet",
                        required=True,
                        choices=['AIS', 'GIS'],
                        help="ice sheet domain: 'AIS' (Antarctica) or 'GIS' (Greenland)")
    parser.add_argument("--authors", dest="authors",
                        required=False, default=DEFAULT_AUTHORS,
                        help=f"author string for output file metadata "
                             f"(default: '{DEFAULT_AUTHORS}')")
    parser.add_argument("--group", dest="group",
                        required=False, default=DEFAULT_GROUP,
                        help=f"group/institution string for output file metadata "
                             f"(default: '{DEFAULT_GROUP}')")
    parser.add_argument("--group_nickname", dest="group_nickname",
                        required=False, default=DEFAULT_GROUP_NICKNAME,
                        help=f"short group nickname used in output filenames "
                             f"(default: '{DEFAULT_GROUP_NICKNAME}')")
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
    if not args.input_file_state is None or not args.input_file_flux is None:
        # Check grid file and res if 2d variables are to be processed
        check_ismip7_grid_file(args.ismip7_grid_file, args.res_ismip7_grid)
        # Check mapping method and either reuse an existing map or create a new one.
        # Note: the function 'building_mapping_file' requires the mpas mesh tool
        # script 'create_SCRIP_file_from_planar_rectangular_grid.py'
        method_remap = args.method_remap
        if args.reuse_mapping_file is not None:
            if not os.path.exists(args.reuse_mapping_file):
                raise FileNotFoundError(f"Mapping file to reuse not found: "
                                        f"{args.reuse_mapping_file}")
            print(f"Reusing existing mapping file: {args.reuse_mapping_file}")
            mapping_file = args.reuse_mapping_file
        else:
            if args.input_file_mesh is None:
                raise ValueError("--input_mesh is required when creating a new "
                                 "mapping file.")

            created_at = datetime.now().strftime("%Y%m%dT%H%M%S")
            mapping_file = f"mapping_mali_to_ismip7.{method_remap}.{created_at}.nc"

            print(f"Creating new mapping file."
                  f"Mapping method used: {method_remap}")

            build_mapping_file(args.input_file_mesh, mapping_file,
                               args.res_ismip7_grid, args.ismip7_grid_file,
                               method_remap)

    print("---Processing remapping file complete---\n")

    # define the path to which the output (processed) files will be saved
    output_path = args.output_path
    print(f"Using output path: {output_path}")
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    # process 1D variables
    if args.global_stats_pattern is None:
        print("--- No global stats pattern provided; skipping 1D variable processing.")
    else:
        print("\n---Processing global stats file(s)---")
        global_stats_files = sorted(glob.glob(args.global_stats_pattern))
        check_global_stats_files(global_stats_files)
        generate_output_1d_vars(global_stats_files, output_path, metadata)
        print("---Processing global stats file(s) complete---\n")

    # process 2d state variables
    if args.input_file_state is None:
        print("--- MALI state file is not provided, thus it will not be processed.")
    else:
        print("\n---Processing state file---")
        # state variables processing part
        # process (add and rename) state vars as requested by the ISMIP7 protocol
        print("Calculating needed state file adjustments.")
        tmp_file = "tmp_state.nc"
        process_state_vars(args.input_file_state, tmp_file)

        # remap data from the MALI unstructured mesh to the ISMIP7 polarstereo grid
        processed_and_remapped_file_state = f'processed_and_remapped_' \
                               f'{os.path.basename(args.input_file_state)}'

        print("Remapping state file.")
        command = ["ncremap",
                   "-i", tmp_file,
                   "-o", processed_and_remapped_file_state,
                   "-m", mapping_file,
                   "-P", "mpas"]
        check_call(command)

        # write out 2D state output files in the ismip7-required format
        print("Writing processed and remapped state fields to ISMIP7 file format.")
        generate_output_2d_state_vars(processed_and_remapped_file_state,
                                      args.ismip7_grid_file,
                                      output_path, metadata)

        os.remove(tmp_file)
        os.remove(processed_and_remapped_file_state)
        print("---Processing state file complete---\n")

    # process 2d flux variables
    if args.input_file_flux is None:
        print("--- MALI flux file is not provided, thus it will not be processed.")
    else:
        print("\n---Processing flux file---")

        # remap data from the MALI unstructured mesh to the ISMIP7 P-S grid
        processed_file_flux = f'processed_' \
                              f'{os.path.basename(args.input_file_flux)}'
        command = ["ncremap",
                   "-i", args.input_file_flux,
                   "-o", processed_file_flux,
                   "-m", mapping_file,
                   "-P", "mpas"]
        check_call(command)

        # write out the output files in the ismip7-required format
        generate_output_2d_flux_vars(processed_file_flux,
                                     args.ismip7_grid_file,
                                     output_path, metadata)

        cleanUp = True
        if cleanUp:
            os.remove(processed_file_flux)
        print("---Processing flux file complete---\n")
    print("---All processing complete---")

if __name__ == "__main__":
    main()
