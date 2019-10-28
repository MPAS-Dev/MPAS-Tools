#!/usr/bin/env python
"""
This script creates coupling files needed to run an MPAS-Ocean and
MPAS-seaice within E3SM.
"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
import os
import shutil
import subprocess
from six.moves.configparser import ConfigParser
import numpy as np


def main():

    config = ConfigParser()
    config.read("create_E3SM_MPAS_coupling_files.ini")

    initial_condition_ocean(config)
    graph_partition_ocean(config)
    initial_condition_seaice(config)
    scrip(config)
    transects_and_regions(config)
    mapping_Gcase(config)
    domain(config)


def initial_condition_ocean(config):
    section = 'initial_condition_ocean'
    if config.get(section, 'enable') == False:
        return
    make_dir(section)
    os.chdir(section)

    # obtain configuration settings
    mesh_name = config.get('main', 'mesh_name')
    coupling_output = config.get('main', 'coupling_file_output')

    # create links
    make_link('../' + mesh_name + '.nc', mesh_name + '.nc')

    # command line execution
    args = ['ncks', '-x', '-v', 'xtime',
            '-O',  # overwrite file if it exists
            mesh_name + '.nc',
            mesh_name + '_no_xtime.nc'
            ]

    run_command(args)

    # create link to output directory
    make_dir('../' + coupling_output + '/ocn/mpas-o/' + mesh_name)
    os.chdir('../' + coupling_output + '/ocn/mpas-o/' + mesh_name)
    make_link(
    '../../../../../' +
    section +
    '/' +
    mesh_name +
    '_no_xtime.nc',
    mesh_name +
     '_no_xtime.nc')
    os.chdir('../../../../..')


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
    mesh_name = config.get('main', 'mesh_name')
    date_string = config.get('main', 'date_string')
    coupling_output = config.get('main', 'coupling_file_output')

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
    make_dir('../' + coupling_output + '/ocn/mpas-o/' + mesh_name)
    os.chdir('../' + coupling_output + '/ocn/mpas-o/' + mesh_name)
    for file in os.listdir('../../../../../' + section):
        if file.startswith('mpas-o.graph.info.' + date_string + '.part.'):
            make_link('../../../../../' + section + '/' + file, './' + file)

    os.chdir('../../../../..')


def initial_condition_seaice(config):
    section = 'initial_condition_seaice'
    if config.get(section, 'enable') == False:
        return
    make_dir(section)
    os.chdir(section)

    # obtain configuration settings
    mesh_name = config.get('main', 'mesh_name')
    coupling_output = config.get('main', 'coupling_file_output')

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
    make_dir('../' + coupling_output + 'ice/mpas-cice/' + mesh_name)
    os.chdir('../' + coupling_output + 'ice/mpas-cice/' + mesh_name)
    make_link(
    '../../../../../' +
    section +
    '/seaice.' +
    mesh_name +
    '.nc',
    'seaice.' +
    mesh_name +
     '.nc')

    os.chdir('../../../../../')


def scrip(config):
    section = 'scrip'
    if config.get(section, 'enable') == False:
        return
    make_dir(section)
    os.chdir(section)

    # obtain configuration settings
    mesh_name = config.get('main', 'mesh_name')
    MPAS_Tools = config.get('main', 'MPAS_Tools')
    date_string = config.get('main', 'date_string')
    mesh_name = config.get('main', 'mesh_name')
    coupling_output = config.get('main', 'coupling_file_output')

    # name scrip files
    scrip_file = (mesh_name + '.nc')
    ocean_scrip_file = ('ocean.' + mesh_name + '.scrip.' + date_string + '.nc')

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
        '-m',
        scrip_file,
        '-s',
        ocean_scrip_file]

    run_command(args)

    # make links to output directories
    make_dir('../' + coupling_output + 'ocn/mpas-o/' + mesh_name)
    os.chdir('../' + coupling_output + 'ocn/mpas-o/' + mesh_name)
    make_link('../../../../../scrip/' + ocean_scrip_file, ocean_scrip_file)
    os.chdir('../../../../..')


def transects_and_regions(config):
    section = 'transects_and_regions'
    if config.get(section, 'enable') == False:
        return
    make_dir(section)
    os.chdir(section)

    # obtain configuration settings
    MPAS_Tools = config.get('main', 'MPAS_Tools')
    mesh_name = config.get('main', 'mesh_name')
    date_string = config.get('main', 'date_string')
    coupling_output = config.get('main', 'coupling_file_output')

    maskfile = 'masks_SingleRegionAtlanticWTransportTransects.' + mesh_name + '.nc'

    # make links
    make_link(
        '/usr/projects/climate/mpeterse/repos/tools/compiled_mesh_conversion_tools/mesh_tools/mesh_conversion_tools/bin/gr/gnu/MpasMaskCreator.x',
        'MpasMaskCreator.x')
    make_link('../' + mesh_name + '.nc', mesh_name + '.nc')
    make_link(
        '/usr/projects/climate/mpeterse/analysis_input_files/geojson_files/SingleRegionAtlanticWTransportTransects.geojson',
        'SingleRegionAtlanticWTransportTransects.geojson')

    argv = ['./MpasMaskCreator.x', mesh_name + '.nc', maskfile,
            '-f', 'SingleRegionAtlanticWTransportTransects.geojson']
    run_command(argv)

    # make links in output directory
    make_dir('../' + coupling_output + 'ocn/mpas-o/' + mesh_name)
    os.chdir('../' + coupling_output + 'ocn/mpas-o/' + mesh_name)
    make_link('../../../../../' + section + '/masks_SingleRegionAtlanticWTransportTransects.' + mesh_name + '.nc',
             'masks_SingleRegionAtlanticWTransportTransects.' + mesh_name + '.nc')

    os.chdir('../../../../..')


def mapping_Gcase(config):
    section = 'mapping_Gcase'
    if config.get(section, 'enable') == False:
        return
    make_dir(section)
    os.chdir(section)

    # obtain configuration settings
    MPAS_Tools = config.get('main', 'MPAS_Tools')
    mesh_name = config.get('main', 'mesh_name')
    date_string = config.get('main', 'date_string')
    E3SM_input = config.get('main', 'E3SM_input')
    coupling_output = config.get('main', 'coupling_file_output')

    # make links
    make_link(
    '../scrip/ocean.' +
    mesh_name +
    '.scrip.' +
    date_string +
    '.nc',
    'ocean.' +
    mesh_name +
    '.scrip.' +
    date_string +
     '.nc')
    make_link(E3SM_input + 'share/scripgrids/T62_040121.nc', 'T62_040121.nc')

    # execute commands
    try:
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--source',
       'ocean.' + mesh_name + '.scrip.' + date_string + '.nc',
       '--destination', 'T62_040121.nc', '--method', 'conserve', '--weight',
       'map_' + mesh_name + '_TO_T62_040121_aave.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--source',
       'ocean.' + mesh_name + '.scrip.' + date_string + '.nc',
       '--destination', 'T62_040121.nc', '--method', 'bilinear', '--weight',
       'map_' + mesh_name + '_TO_T62_040121_blin.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--source',
       'ocean.' + mesh_name + '.scrip.' + date_string + '.nc',
       '--destination', 'T62_040121.nc', '--method', 'patch', '--weight',
       'map_' + mesh_name + '_TO_T62_040121_patc.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--destination',
       'ocean.' + mesh_name + '.scrip.' + date_string + '.nc',
       '--source', 'T62_040121.nc', '--method', 'conserve', '--weight',
       'map_T62_040121_TO_' + mesh_name + '_aave.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--destination',
       'ocean.' + mesh_name + '.scrip.' + date_string + '.nc',
       '--source', 'T62_040121.nc', '--method', 'bilinear', '--weight',
       'map_T62_040121_TO_' + mesh_name + '_blin.' + date_string + '.nc',
       '--ignore_unmapped']
       run_command(argv)
       
       argv = ['mpirun', '-n', '1', 'ESMF_RegridWeightGen', '--destination',
       'ocean.' + mesh_name + '.scrip.' + date_string + '.nc',
       '--source', 'T62_040121.nc', '--method', 'patch', '--weight',
       'map_T62_040121_TO_' + mesh_name + '_patc.' + date_string + '.nc',
       '--ignore_unmapped'] 
       run_command(argv)

    except OSError: 
        print('mapping_Gcase must be run on one compute node')
   
    # make links in output directory
    make_dir('../' + coupling_output + 'cpl/cpl6')
    os.chdir('../' + coupling_output + 'cpl/cpl6')
    for file in os.listdir('../../../../mapping_Gcase'):
        if file.startswith('map_'):
            make_link('../../../../mapping_Gcase/' + file, './' + file)

    os.chdir('../../../..')

def domain(config):
    section = 'domain'
    if config.get(section, 'enable') == False:
        return
    make_dir(section)
    os.chdir(section)

    # obtain configuration settings
    date_string = config.get('main', 'date_string')
    mesh_name = config.get('main', 'mesh_name')
    E3SM_data = config.get('main', 'E3SM_data_repo')
    coupling_output = config.get('main', 'coupling_file_output')

    # make links
    make_link(E3SM_data + 'compiled_cime_tools/cime/tools/mapping/gen_domain_files/src/gen_domain', 'gen_domain')
    make_link('../mapping_Gcase/map_' + mesh_name +'_TO_T62_040121_aave.' + date_string +'.nc',
              'map_' + mesh_name +'_TO_T62_040121_aave.' + date_string +'.nc')
    
# execute commands
    argv = ['./gen_domain', '-m', 'map_' + mesh_name + '_TO_T62_040121_aave.' + date_string + '.nc', '-o', mesh_name, '-l', 'T62']
    run_command(argv)

    # make links in output directories
    make_dir('../' + coupling_output + 'share/domains')
    os.chdir('../' + coupling_output + '/share/domains')
    make_link('../../../../' + section + '/domain.lnd.T62_' + mesh_name + '.' + date_string + '.nc',
              'domain.lnd.T62_' + mesh_name + '.' + date_string + '.nc' )
    make_link('../../../../' + section + '/domain.ocn.T62_' + mesh_name + '.' + date_string + '.nc',
              'domain.ocn.T62_' + mesh_name + '.' + date_string + '.nc')
    make_link('../../../../' + section + '/domain.ocn.' + mesh_name + '.' + date_string + '.nc',
              'domain.ocn' + mesh_name + '.' + date_string + '.nc')
    
    os.chdir('../../../..')

def make_dir(dirName):
    try:
        write_command_history(dirName)
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
        write_command_history('  ' + ' '.join(args))
        with open('log.out', 'a') as outstream:
            outstream.write('Command: ' + ' '.join(args) + '\n')
            subprocess.check_call(args, stdout=outstream, stderr=outstream)
            outstream.write('\n')
    except OSError:
        pass


if __name__ == '__main__':
    # If called as a primary module, run main
    main()
