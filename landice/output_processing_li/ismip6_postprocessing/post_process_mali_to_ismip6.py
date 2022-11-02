#!/usr/bin/env python

"""
This script processes MALI simulation outputs (both state and flux)
in the required format by the ISMIP6 experimental protocol.
The input state files (i.e., output files from MALI) need to have been
concatenated to have yearly data, which can be done using 'ncrcat' command
before using this script.
"""

import argparse
from subprocess import check_call
import os
import shutil
from create_mapfile_mali_to_ismip6 import build_mapping_file
from process_state_variables import generate_output_2d_state_vars, \
     process_state_vars, generate_output_1d_vars
from process_flux_variables import generate_output_2d_flux_vars, \
     do_time_avg_flux_vars, clean_flux_fields_before_time_averaging


def main():
    parser = argparse.ArgumentParser(
                        description='process MALI outputs for the ISMIP6'
                                    'submission')
    parser.add_argument("-e", "--exp_name", dest="exp",
                        required=True,
                        help="ISMIP6 experiment name (e.g., exp05")
    parser.add_argument("-i_state", "--input_state", dest="input_file_state",
                        required=False, help="mpas output state variables")
    parser.add_argument("-i_flux", "--input_flux", dest="input_file_flux",
                        required=False, help="mpas output flux variables")
    parser.add_argument("-i_init", "--input_grid", dest="input_file_grid",
                        required=False, help="MALI file with mesh information"
                                            "Only used for creating mapping file,"
                                            "So not required if mapping file already exists"
                                            "Can be a state or flux output file.")
    parser.add_argument("-g", "--global_stats_file", dest="global_stats_file",
                        required=False, help="globalStats.nc file")
    parser.add_argument("-p", "--output_path", dest="output_path",
                        required=False,
                        help="path to which the final output files"
                             " will be saved")
    parser.add_argument("--mali_mesh_name", dest="mali_mesh_name",
                        required=True,
                        help="mali mesh name (e.g., AIS_8to30km")
    parser.add_argument("--mapping_file", dest="mapping_file",
                        required=False,
                        help="mapping file name from MALI mesh to ISMIP6 grid")
    parser.add_argument("--ismip6_grid_file", dest="ismip6_grid_file",
                        required=True,
                        help="Input ismip6 mesh file.")
    parser.add_argument("--method", dest="method_remap", default="conserve",
                        required=False,
                        help="mapping method. Default='conserve'")

    args = parser.parse_args()

    print("\n---Processing remapping file---")
    # Only do remapping steps if we have 2d files to process
    if not args.input_file_state is None or not args.input_file_flux is None:
        # check the mapping method and existence of the mapping file
        # Note: the function 'building_mapping_file' requires the mpas mesh tool
        # script 'create_SCRIP_file_from_planar_rectangular_grid.py'
        if os.path.exists(args.mapping_file):
            print(f"Mapping file exists.")
            mapping_file = args.mapping_file
        else:
            if args.method_remap is None:
                method_remap = "conserve"
            else:
                method_remap = args.method_remap

            mapping_file = f"map_{args.mali_mesh_name}_to_ismip6_8km_" \
                           f"{method_remap}.nc"

            print(f"Creating new mapping file."
                  f"Mapping method used: {method_remap}")
            build_mapping_file(args.input_file_grid, mapping_file,
                               args.ismip6_grid_file, method_remap)

    print("---Processing remapping file complete---\n")

    # define the path to which the output (processed) files will be saved
    if args.output_path is None:
        output_path = os.getcwd()
    else:
        output_path = args.output_path
    print(f"Using output path: {output_path}")
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    if args.input_file_state is None:
        print("--- MALI state file is not provided, thus it will not be processed.")
    else:
        print("\n---Processing state file---")
        # state variables processing part
        # process (add and rename) state vars as requested by the ISMIP6 protocol
        print("Calculating needed state file adjustments.")
        tmp_file = "tmp_state.nc"
        process_state_vars(args.input_file_state, tmp_file)

        # remap data from the MALI unstructured mesh to the ISMIP6 polarstereo grid
        processed_and_remapped_file_state = f'processed_and_remapped_' \
                               f'{os.path.basename(args.input_file_state)}'

        print("Remapping state file.")
        command = ["ncremap",
                   "-i", tmp_file,
                   "-o", processed_and_remapped_file_state,
                   "-m", mapping_file,
                   "-P", "mpas"]
        check_call(command)

        # write out 2D state output files in the ismip6-required format
        print("Writing processed and remapped state fields to ISMIP6 file format.")
        generate_output_2d_state_vars(processed_and_remapped_file_state, args.ismip6_grid_file,
                                      args.exp, output_path)

        os.remove(tmp_file)
        os.remove(processed_and_remapped_file_state)
        print("---Processing state file complete---\n")

    # write out 1D output files for both state and flux variables
    if args.global_stats_file is None:
        print("--- MALI global stats file is not provided, thus it will not be processed.")
    else:
        print("\n---Processing global stats file---")
        generate_output_1d_vars(args.global_stats_file, args.exp,
                                output_path)
        print("---Processing global stats file complete---\n")

    # process the flux variables if flux output file is given
    if args.input_file_flux is None:
        print("--- MALI flux file is not provided, thus it will not be processed.")
    else:
        print("\n---Processing flux file---")

        print("Translate calving flux")
        tmp_file_translate = "flux_translated.nc"
        clean_flux_fields_before_time_averaging(args.input_file_flux, tmp_file_translate)
        # take time (yearly) average for the flux variables
        tmp_file1 = "flux_time_avg.nc"
        do_time_avg_flux_vars(tmp_file_translate, tmp_file1)

        # remap data from the MALI unstructured mesh to the ISMIP6 P-S grid
        processed_file_flux = f'processed_' \
                              f'{os.path.basename(args.input_file_flux)}'
        command = ["ncremap",
                   "-i", tmp_file1,
                   "-o", processed_file_flux,
                   "-m", mapping_file,
                   "-P", "mpas"]
        check_call(command)

        # write out the output files in the ismip6-required format
        generate_output_2d_flux_vars(processed_file_flux,
                                     args.ismip6_grid_file,
                                     args.exp, output_path)

        cleanUp = True
        if cleanUp:
            os.remove(tmp_file_translate)
            os.remove(tmp_file1)
            os.remove(processed_file_flux)
        print("---Processing flux file complete---\n")
    print("---All processing complete---")

if __name__ == "__main__":
    main()
