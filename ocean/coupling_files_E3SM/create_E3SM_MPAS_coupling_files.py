#!/usr/bin/env python
"""
This script creates coupling files needed to run an MPAS-Ocean and
MPAS-seaice within E3SM.

Load the lastest e3sm-unified conda package.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
import os
import shutil
import subprocess
from six.moves.configparser import ConfigParser
import numpy as np
import glob


def main():

    # obtain configuration settings
    config = ConfigParser()
    config.read("create_E3SM_MPAS_coupling_files.ini")
    mesh_name = config.get('main', 'mesh_name')

    function_list = ['initial_condition_ocean',
                     'graph_partition_ocean',
                     'initial_condition_seaice',
                     'scrip',
                     'transects_and_regions',
                     'mapping_Gcase',
                     'domain'];

    # create input_data directories
    make_dir('assembled_files_for_upload/input_data/ocn/mpas-o/' + mesh_name)
    make_dir('assembled_files_for_upload/input_data/ice/mpas-cice/' + mesh_name)
    make_dir('assembled_files_for_upload/input_data/cpl/cpl6')
    make_dir('assembled_files_for_upload/input_data/share/domains')

    for function_name in function_list:
        print("****** " + function_name + " ******")

        if config.get(function_name, 'enable') == False:
            print("Disabled in .ini file")
            return
        currentDir = os.getcwd()
        make_dir(function_name)
        os.chdir(function_name)

        try:
            globals()[function_name](function_name, config)
            print('SUCCESS: ' + function_name + ' completed.')
        except:
            print('WARNING: ' + function_name + ' failed.')
        print(" ")
        os.chdir(currentDir)


def initial_condition_ocean(function_name,config):

    # obtain configuration settings
    mesh_name = config.get('main', 'mesh_name')

    # create links
    make_link('../' + mesh_name + '.nc', mesh_name + '.nc')

    # command line execution
    args = ['ncks', '-x', '-v', 'xtime', '-O',
            mesh_name + '.nc',
            mesh_name + '_no_xtime.nc'
            ]
    run_command(args)

    # create link to output directory
    os.chdir('../assembled_files_for_upload/input_data/ocn/mpas-o/' + mesh_name)
    make_link('../../../../../' + function_name + '/' + mesh_name + '_no_xtime.nc',
        mesh_name + '_no_xtime.nc')


def graph_partition_ocean(function_name,config):

    # obtain configuration settings
    min_graph_size = config.getint(function_name, 'min_graph_size')
    max_graph_size = config.getint(function_name, 'max_graph_size')
    graph_file = config.get(function_name, 'graph_file')
    mesh_name = config.get('main', 'mesh_name')
    date_string = config.get('main', 'date_string')

    # create links
    make_link('../' + graph_file, 'mpas-o.graph.info.' + date_string)

    # command line execution
    n_power2 = 2**np.arange(3, 10 + 1)
    n_multiples12 = 12 * np.arange(1, 8 + 1)

    n = n_power2
    for power10 in range(3):
        n = np.concatenate([n, 10**power10 * n_multiples12])

    for j in range(len(n)):
        if n[j] >= min_graph_size and n[j] <= max_graph_size:
            args = ['gpmetis', 'mpas-o.graph.info.' + date_string, str(n[j])]
            run_command(args)

    # create link to output directory
    files = glob.glob('mpas-o.graph.info.*')
    os.chdir('../assembled_files_for_upload/input_data/ocn/mpas-o/' + mesh_name)
    for file in files:
        make_link('../../../../../graph_partition_ocean/' + file, './' + file)


def initial_condition_seaice(function_name,config):

    # obtain configuration settings
    mesh_name = config.get('main', 'mesh_name')

    # create links
    make_link('../' + mesh_name + '.nc', mesh_name + '.nc')

    # command line execution
    args = ['ncks', '-x', '-v',
            'bottomDepth,refBottomDepth,restingThickness, temperature,salinity,\
    temperatureSurfaceValue,salinitySurfaceValue,surfaceVelocityZonal,\
    surfaceVelocityMeridional,SSHGradientZonal,SSHGradientMeridional,\
    vertNonLocalFluxTemp,normalVelocity,layerThickness,normalBarotropicVelocity,\
    vertCoordMovementWeights,boundaryLayerDepth,seaIcePressure,\
    atmosphericPressure,filteredSSHGradientZonal,filteredSSHGradientMeridional',
            '-O',  # Overwrite existing file
            mesh_name + '.nc',
            'seaice.' + mesh_name + '.nc'
            ]

    run_command(args)

    # make links to output directory
    os.chdir('../assembled_files_for_upload/input_data/ice/mpas-cice/' + mesh_name)
    make_link( '../../../../../' + function_name + '/seaice.' + mesh_name + '.nc',
               'seaice.' + mesh_name + '.nc')


def scrip(function_name,config):

    # obtain configuration settings
    mesh_name = config.get('main', 'mesh_name')
    MPAS_Tools = config.get('main', 'MPAS_Tools')
    date_string = config.get('main', 'date_string')
    mesh_name = config.get('main', 'mesh_name')

    # name scrip files
    scrip_file = ('ocean.' + mesh_name + '.scrip.' + date_string + '.nc')

    # create links
    make_link('../' + mesh_name + '.nc', mesh_name + '.nc')
    make_link(
        MPAS_Tools +
        '/mesh_tools/create_SCRIP_files/create_SCRIP_file_from_MPAS_mesh.py',
        'create_SCRIP_file_from_MPAS_mesh.py')

    # command line execution
    args = [
        'python',
        'create_SCRIP_file_from_MPAS_mesh.py',
        '-m', mesh_name + '.nc',
        '-s', scrip_file]

    run_command(args)

    # make links to output directories
    os.chdir('../assembled_files_for_upload/input_data/ocn/mpas-o/' + mesh_name)
    make_link('../../../../../scrip/' + scrip_file, scrip_file)


def transects_and_regions(function_name,config):

    # obtain configuration settings
    mesh_name = config.get('main', 'mesh_name')
    date_string = config.get('main', 'date_string')
    transect_region_geojson = config.get('transects_and_regions', 'transect_region_geojson')

    maskfile = 'masks_SingleRegionAtlanticWTransportTransects.' + mesh_name + '.nc'

    # make links
    make_link('../' + mesh_name + '.nc', mesh_name + '.nc')
    make_link(transect_region_geojson,
        'SingleRegionAtlanticWTransportTransects.geojson')

    # check: how does this call MpasMaskCreator?
    argv = ['MpasMaskCreator.x', mesh_name + '.nc', maskfile,
            '-f', 'SingleRegionAtlanticWTransportTransects.geojson']
    run_command(argv)

    # make links in output directory
    os.chdir('../assembled_files_for_upload/input_data/ocn/mpas-o/' + mesh_name)
    make_link('../../../../../' + function_name + '/masks_SingleRegionAtlanticWTransportTransects.' + mesh_name + '.nc',
             'masks_SingleRegionAtlanticWTransportTransects.' + mesh_name + '.nc')


def mapping_Gcase(function_name,config):

    # obtain configuration settings
    mesh_name = config.get('main', 'mesh_name')
    date_string = config.get('main', 'date_string')
    E3SM_input_data = config.get('main', 'E3SM_input_data')

    # make links
    scrip_file = 'ocean.' + mesh_name + '.scrip.' + date_string + '.nc'
    make_link('../scrip/ocean.' + mesh_name + '.scrip.' + date_string + '.nc',
              scrip_file)
    make_link(E3SM_input_data + 'share/scripgrids/T62_040121.nc', 'T62_040121.nc')

    # execute commands
    try:
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--source',
       scrip_file,
       '--destination', 'T62_040121.nc', '--method', 'conserve', '--weight',
       'map_' + mesh_name + '_TO_T62_040121_aave.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--source',
       scrip_file,
       '--destination', 'T62_040121.nc', '--method', 'bilinear', '--weight',
       'map_' + mesh_name + '_TO_T62_040121_blin.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--source',
       scrip_file,
       '--destination', 'T62_040121.nc', '--method', 'patch', '--weight',
       'map_' + mesh_name + '_TO_T62_040121_patc.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--destination',
       scrip_file,
       '--source', 'T62_040121.nc', '--method', 'conserve', '--weight',
       'map_T62_040121_TO_' + mesh_name + '_aave.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--destination',
       scrip_file,
       '--source', 'T62_040121.nc', '--method', 'bilinear', '--weight',
       'map_T62_040121_TO_' + mesh_name + '_blin.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--destination',
       scrip_file,
       '--source', 'T62_040121.nc', '--method', 'patch', '--weight',
       'map_T62_040121_TO_' + mesh_name + '_patc.' + date_string + '.nc',
       '--ignore_unmapped'] 
       run_command(argv)

    except OSError: 
        print('mapping_Gcase must be run on one compute node')
   
    # make links in output directory
    files = glob.glob('map_*')
    os.chdir('../assembled_files_for_upload/input_data/cpl/cpl6')
    for file in files:
        make_link('../../../../mapping_Gcase/' + file, './' + file)


def domain(function_name,config):

    # obtain configuration settings
    date_string = config.get('main', 'date_string')
    mesh_name = config.get('main', 'mesh_name')
    gen_domain = config.get('domain', 'gen_domain')

    # make links
    make_link(gen_domain, 'gen_domain')
    mapping_file = 'map_' + mesh_name + '_TO_T62_040121_aave.' + date_string +'.nc'
    make_link('../mapping_Gcase/' + mapping_file, mapping_file)
    
    # execute commands
    argv = ['./gen_domain', '-m', mapping_file, '-o', mesh_name, '-l', 'T62']
    run_command(argv)

    # make links in output directories
    files = glob.glob('domain*.nc')
    os.chdir('../assembled_files_for_upload/input_data/share/domains')
    for file in files:
        make_link('../../../../domain/' + file, './' + file)
    

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


def write_command_history(text):
    try:
        print(text)
        with open('command_history', 'a') as outstream:
            outstream.write(text + '\n')
    except OSError:
        pass


def run_command(args):
    try:
        write_command_history(' '.join(args))
        with open('log.out', 'a') as outstream:
            outstream.write('Command: ' + ' '.join(args) + '\n')
            subprocess.check_call(args, stdout=outstream, stderr=outstream)
            outstream.write('\n')
    except OSError:
        pass


if __name__ == '__main__':
    # If called as a primary module, run main
    main()
