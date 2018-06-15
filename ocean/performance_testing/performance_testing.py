#!/usr/bin/env python

"""
 Name: performance_testing.py
 Author: Divya Jaganathan
 Date: 6 July, 2018

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

os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OMP_PLACES'] = 'threads'

timenow = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

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

cells_fr = open("graph.info", 'r')
cells = str(cells_fr.readline().split(" ")[0])

nprocs_max = int(sys.argv[1])
nprocs_min = int(sys.argv[2])
cpu_type = sys.argv[3]
cores_per_node = float(sys.argv[4])

plane_size = str(int(cores_per_node))

niter = int(np.log2(nprocs_max)) - int(np.log2(nprocs_min)) + 1
nsamples_per_procnum = 5

time = np.zeros(shape=(1, niter))
procs = np.zeros(shape=(1, niter))
SYPD = np.zeros(shape=(1, niter))

i = nprocs_max
j = niter

writefilename = "data_" + cpu_type + "_" + \
    str(nprocs_max) + "_" + timenow + ".txt"
fw = open(writefilename, 'a+')
fw.write(
    'Time: %s \nMachine: %s\nNo. of Cells: %s\nRun duration: %s\nRun time in sec: %d\nFormat: #Procs|Sample Runs|Average|SYPD \n' %
    (timenow, cpu_type, cells, timeparts, simulated_time_in_sec))

while i >= nprocs_min:

    local_N = int(np.ceil(i / cores_per_node))
    sample = nsamples_per_procnum
    foldername = "perf_p" + str(i) + "_gr_openmpi"
    subprocess.check_call(['rm', '-rf', foldername])
    subprocess.check_call(['mkdir', foldername])
    fw.write('%s \t' % i)
    sum = 0

    # Generate the log and graph files
    subprocess.check_call(['./metis', 'graph.info', str(i)])
    print "metis" + str(i) + "completed"

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

# plotting ..

subprocess.check_call(['mkdir', '-p', 'data_figures'])

perfect = SYPD[0][0] / procs[0][0] * procs
plt.loglog(procs[0], SYPD[0], '-or', label='grizzly')
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

# End Version
