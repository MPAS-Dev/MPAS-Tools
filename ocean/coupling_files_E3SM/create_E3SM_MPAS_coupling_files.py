#!/usr/bin/env python
"""
This script creates coupling files needed to run an MPAS-Ocean and
MPAS-seaice within E3SM.
"""
import os
import shutil
import subprocess
import ConfigParser
import numpy as np


def main():

    config = ConfigParser.ConfigParser()
    config.read("create_E3SM_MPAS_coupling_files.ini")

    initial_condition_ocean(config)
    graph_partition_ocean(config)


def initial_condition_ocean(config):
    section = 'initial_condition_ocean'
    if config.get(section, 'enable') == False:
        return
    make_dir(section)
    os.chdir(section)

    # obtain configuration settings
    input_file = config.get('main', 'input_file')

    # create links
    make_link('../' + input_file, input_file)

    # command line execution
    args = ['ncks', '-x', '-v', 'xtime',
            input_file,
            input_file + '_no_xtime.nc'
            ]

    # run_command(args)
    os.chdir('..')


def graph_partition_ocean(config):
    section = 'graph_partition_ocean'
    if config.get(section, 'enable') == False:
        return
    make_dir(section)
    os.chdir(section)

    # obtain configuration settings
    min_graph_size = config.getint(section, 'min_graph_size')
    max_graph_size = config.getint(section, 'max_graph_size')
    graph_file = config.get(section, 'graph_file')
    gpmetis_executable = config.get(section, 'gpmetis_executable')
    date_string = config.get('main', 'date_string')

    # create links
    make_link('../' + graph_file, 'mpas-o.graph.info.' + date_string)
    make_link(gpmetis_executable, 'gpmetis')

    # command line execution
    n_power2 = 2**np.arange(3, 10 + 1)
    n_multiples12 = 12 * np.arange(1, 8 + 1)

    n = n_power2
    for power10 in range(3):
        n = np.concatenate([n, 10**power10 * n_multiples12])

    for j in range(len(n)):
        if n[j] >= min_graph_size and n[j] <= max_graph_size:
            args = ['./gpmetis', 'mpas-o.graph.info.' + date_string, str(n[j])]
            run_command(args)

    os.chdir('..')


def make_dir(dirName):
    try:
        os.makedirs(dirName)
    except OSError:
        pass


def make_link(source, linkName):
    try:
        if os.path.exists(linkName):
            os.remove(linkName)
        os.symlink(source, linkName)
    except OSError:
        pass


def run_command(args):
    try:
        print(' '.join(args))
        with open('log.out', 'a') as outstream:
            outstream.write('Command: ' + ' '.join(args) + '\n')
            subprocess.check_call(args, stdout=outstream, stderr=outstream)
            outstream.write('\n')
    except OSError:
        pass


if __name__ == '__main__':
    # If called as a primary module, run main
    main()
