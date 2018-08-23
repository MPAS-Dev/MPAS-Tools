#!/usr/bin/env python

"""
 Name: performance_testing.py
 Author: Divya Jaganathan
 Date: 6 July, 2018

This script is automatically called by call_to_performance_testing.py to run a batch job to get performance plots and data

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

# To obtain the run duration used in calculating SYPD

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

# To store the details on number of cells (~ resolution)

cells_fr = open("graph.info", 'r')
cells = str(cells_fr.readline().split(" ")[0])
cells_fr.close()

nprocs_max = int(sys.argv[1])
nprocs_min = int(sys.argv[2])
cpu_type = sys.argv[3]
cores_per_node = float(sys.argv[4])
x = float(sys.argv[5])

# plane_size is used to define the plane_distribution flag in srun
plane_size = str(int(cores_per_node))

# Performance data evaluation begins here -

niter = int(np.log2(nprocs_max)) - int(np.log2(nprocs_min)) + 1
nsamples_per_procnum = 5

time = np.zeros(shape=(1, niter))
procs = np.zeros(shape=(1, niter))
SYPD = np.zeros(shape=(1, niter))

#i = nprocs_max
#j = niter

writefilename = "data_" + cpu_type + "_" + \
    str(nprocs_max) + "_" + timenow + ".txt"
fw = open(writefilename, 'a+')
fw.write(
    'Time: %s \nMachine: %s\nNo. of Cells: %s\nRun duration: %s\nRun time in sec: %d\nFormat: #Procs|Sample Runs|Average|SYPD \n' %
    (timenow, cpu_type, cells, timeparts, simulated_time_in_sec))


while x >= 0:

    graph_call = "python generate_graph.info_with_wgts.py -d init.nc -g graph.info -x %s" % x
    g_args = shlex.split(graph_call)
    print "running", ''.join(g_args)
    subprocess.check_call(g_args)

    foldername_wgt = "weight" + str(x)
    subprocess.check_call(['mkdir', foldername_wgt])
    graph_filename = "graph.info_with_wgts_" + str(x)
    i = nprocs_max
    j = niter

    while i >= nprocs_min:

        local_N = int(np.ceil(i / cores_per_node))
        sample = nsamples_per_procnum
        foldername = foldername_wgt + "/perf_p" + str(i) + "_gr_openmpi"
        subprocess.check_call(['mkdir', '-p', foldername])
        fw.write('%s \t' % i)
        sum = 0
        subprocess.check_call(['./metis', graph_filename, str(i)])
        print "metis" + str(i) + "completed"
        graph_part_name = graph_filename + ".part." + str(i)
        to_name = "graph.info.part." + str(i)
        subprocess.check_call(['mv', graph_part_name, to_name])

        while sample >= 1:
            args = ['srun',
                    '-N',
                    str(local_N),
                    '-n',
                    str(i),
                    '--cpu_bind=verbose,core',
                    '--distribution=plane=%s' % plane_size,
                    './ocean_model']
            print "running", ''.join(args)
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

            fname = "log_p" + str(i) + "_s" + str(sample)
            filepath = foldername + "/" + fname
            sample = sample - 1
            subprocess.check_call(['mv', 'log.ocean.0000.out', filepath])

        average = sum / nsamples_per_procnum
        time[0][j - 1] = average
        procs[0][j - 1] = i
        SYPD[0][j - 1] = simulated_time_in_sec / (365 * average)
        fw.write('%s \t %s\n' % (str(average), str(SYPD[0][j - 1])))
        i = i / 2
        j = j - 1
        subprocess.check_call(['mv', to_name, foldername])
    x = x - 0.5

# plotting ..

subprocess.check_call(['mkdir', '-p', 'data_figures'])

perfect = SYPD[0][0] / procs[0][0] * procs
plt.loglog(procs[0], SYPD[0], '-or', label=str(x))
plt.loglog(procs[0], perfect[0], '--k', label='perfect scaling')
plt.title(r'MPAS-Ocean Performance Curve (Broadwell 36-cores No HT)')
plt.xlabel('Number of MPI ranks')
plt.ylabel('Simulated Years Per Day (SYPD)')
plt.legend(loc='upper left')
plt.grid()
plt.xlim((1, nprocs_max * 2))
plt.tight_layout()
figurenamepath = "data_figures/fig_" + cpu_type + \
    str(nprocs_max) + "_" + timenow + ".png"
plt.savefig(figurenamepath)
subprocess.check_call(['mv', writefilename, 'data_figures'])

fr.close()
fw.close()

# End
