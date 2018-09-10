#!/usr/bin/env python

"""
 Name: performance_test.py
 Author: Divya Jaganathan
 Date: 6 July, 2018

This script is automatically called by submit_performance_test_to_queue.py to run a batch job to get performance plots and data

This script can also be used for an interactive job submission using the following command format:

command format (to run an interactive job) :
./performance_test.py <Maximum Tasks> <Minimum Tasks> <machine_long_name> <cores_per_node> <resolution_name>

Access files required to run this script:
 1. namelist.ocean
 2. graph.info
 3. metis file (rename gpmetis to metis or vice-versa in this script when creating a soft link)
 4. ocean_model (executable file)

NOTE: When running a large number of tasks (>10k), check the name of log.ocean.0000.out file generated - no. of zeros in the file name changes
"""

import subprocess
import numpy as np
import re
import sys
import datetime
from time import strftime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import shlex

# Setting OMP variables for NO multithreading

os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OMP_PLACES'] = 'threads'

timenow = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

# Read namelist to obtain the run duration used in calculating SYPD
time_fr = open("namelist.ocean", 'r')

for line in time_fr:
    m1 = re.search("config_run_duration", line)
    if m1:
        parts = line.split("=", 1)[1]
        subparts = re.split(':|\'|_', parts)
        timeparts = subparts[2:5]
        dateparts = re.split('-', subparts[1])
        if len(dateparts) == 1:
            simulated_time_in_sec = int(
                timeparts[0]) * 60 * 60 + int(timeparts[1]) * 60 + int(timeparts[2])
        else:
            simulated_time_in_sec = int(dateparts[2]) * 24 * 60 * 60 + int(
                timeparts[0]) * 60 * 60 + int(timeparts[1]) * 60 + int(timeparts[2])

time_fr.close()

subprocess.check_call(['mkdir', '-p', 'data_performance'])
subprocess.check_call(['mkdir', '-p', 'figures_performance'])

# Store the details on number of cells for this resolution

cells_fr = open("graph.info", 'r')
cells = str(cells_fr.readline().split(" ")[0])
cells_fr.close()

nprocs_max = int(sys.argv[1])
nprocs_min = int(sys.argv[2])
machine_long_name = sys.argv[3]
cores_per_node = float(sys.argv[4])
resolution_name = sys.argv[5]

# plane_size is used to define the plane_distribution flag in srun
plane_size = str(int(cores_per_node))

# Performance data evaluation begins here
niter = int(np.log2(nprocs_max)) - int(np.log2(nprocs_min)) + 1
nsamples_per_procnum = 3

time = np.zeros(shape=(1, niter))
procs = np.zeros(shape=(1, niter))
SYPD = np.zeros(shape=(1, niter))

i = nprocs_max
j = niter - 1

writefilename = "data_performance/" + machine_long_name + "_" + \
    str(nprocs_max) + "_" + timenow + ".txt"
fw = open(writefilename, 'a+')
fw.write(
    'Time: %s \nMachine: %s\nResolution: %s\nHorizontal cells: %s\nRun duration: %s\nRun time in sec: %d\nFormat: #Procs|Sample Runs|Average|SYPD \n' %
    (timenow,
     machine_long_name,
     resolution_name,
     cells,
     timeparts,
     simulated_time_in_sec))
fw.flush()

while i >= nprocs_min:

    local_N = int(np.ceil(i / cores_per_node))
    foldername = "perf_p" + str(i) + "_gr_openmpi"
    subprocess.check_call(['rm', '-rf', foldername])
    subprocess.check_call(['mkdir', foldername])
    fw.write('%s \t' % i)
    sum = 0

    # Generate the log and graph files
    subprocess.check_call(['./metis', 'graph.info', str(i)])
    print "metis" + str(i) + "completed"

    for sample in range(nsamples_per_procnum):
        subprocess.check_call(
            ['rm', '-rf', 'log*', 'analysis_members', 'output.nc'])
        args = ['srun',
                '-N',
                str(local_N),
                '-n',
                str(i),
                '--cpu_bind=verbose,core',
                '--distribution=plane=%s' % plane_size,
                './ocean_model']
        print "running", ' '.join(args)
        subprocess.check_call(args)

        # Search for time integration and write to a file
        fr = open("log.ocean.0000.out", 'r')
        for line in fr:
            m = re.search("2  time integration", line)
            if m:
                numbers = line.split("integration", 1)[1]
                first_number = numbers.split()[0]
                fw.write('%s \t' % first_number)
                sum = sum + float(first_number)
        fr.close()
        fname = "log_p" + str(i) + "_s" + str(sample + 1)
        filepath = foldername + "/" + fname
        subprocess.check_call(['mv', 'log.ocean.0000.out', filepath])

    average = sum / nsamples_per_procnum
    time[0][j] = average
    procs[0][j] = i
    SYPD[0][j] = simulated_time_in_sec / (365 * average)
    fw.write('%s \t %s\n' % (str(average), str(SYPD[0][j])))
    fw.flush()
    perfect = SYPD[0][j] / procs[0][j] * procs

    # create plot with data so far
    plt.clf()
    plt.loglog(procs[0][j:], SYPD[0][j:], '-or',
               label=resolution_name + ', ' + machine_long_name)
    plt.loglog(procs[0][j:], perfect[0][j:], '--k', label='perfect scaling')
    plt.title('MPAS-Ocean Performance Curve')
    plt.xlabel('Number of MPI ranks')
    plt.ylabel('Simulated Years Per Day (SYPD)')
    plt.legend(loc='upper left')
    plt.grid(which='major')
    plt.xlim((procs[0][j] / 2.0, nprocs_max * 2.0))
    plt.tight_layout()
    figurenamepath = "figures_performance/" + resolution_name + '_' + \
        machine_long_name + '_' + str(nprocs_max) + "_" + timenow + ".png"
    plt.savefig(figurenamepath)
    i = i / 2
    j = j - 1
fw.close()
