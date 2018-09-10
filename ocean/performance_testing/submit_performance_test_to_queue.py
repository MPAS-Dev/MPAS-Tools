#!/usr/bin/env python

"""
Name: submit_performance_test_to_queue.py
Author: Divya Jaganathan
Date: July 6, 2018

Submits request for a batch job to carry out successive performance runs starting from maximum
number of tasks.  Load modules before calling this script.

command format:
./submit_performance_test_to_queue.py -M <Maximum Tasks> -m <Minimum Tasks> -n <machine_name> -r <resolution_name>

Examples:

On any machine, you can grab tarred run directories here:
https://zenodo.org/record/1252437#.W5FIppNKjUI
add a link to
- metis
- ocean_model executable
- submit_performance_test_to_queue.py (here)
- performance_test.py (here)

On any machine log-in node, all you need is:
   ./submit_performance_test_to_queue.py
This will submit a single job to the queue, and produce the default test of
64 through 2 by powers of 2, and auto-detect your machine.  Load modules before
calling this script, and submission will keep the same modules.

Or, one can specify everything with flags.  This tests 128 to 16 cores by powers of two.
   ./submit_performance_test_to_queue.py -M 128 -m 16 -r EC60to30

On cori, you have to specify cori-knl or cori-haswell, as follows:
   ./submit_performance_test_to_queue.py -M 128 -m 16 -r EC60to30 -n cori-knl

After the job completes, you will find data and auto-generated plots in these directories:
  data_performance
  figures_performance
"""

import subprocess
import argparse
import shlex
import numpy as np
import os

parser = \
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "-M",
    "--max_tasks",
    dest="max_tasks",
    help="Maximum number of tasks, defaults to 64.",
    default=64)
parser.add_argument(
    "-m",
    "--min_tasks",
    dest="min_tasks",
    help="Minimum number of tasks, defaults to 2.",
    default=2)
parser.add_argument(
    "-n",
    "--machine_name",
    dest="machine_name",
    help="This script auto-detects the machine from the node name (e.g. 'gr' for grizzly). Use this flag to override. On cori, enter cori-haswell or cori-knl",
    default=os.uname()[1][0:2])
parser.add_argument(
    "-r",
    "--resolution_name",
    dest="resolution_name",
    help="This label appears on the title of the plot.",
    default="MPAS-O")
args = parser.parse_args()

max_tasks = int(args.max_tasks)
min_tasks = int(args.min_tasks)
machine_name = args.machine_name
resolution_name = args.resolution_name

job_id = "MPASO_perf_P" + str(max_tasks) + args.resolution_name
output_name = "slurm_" + job_id + ".out"

# NODES_REQUIRED to request for resources is calculated assuming no hyperthreads.
# Changes to this can be implemented by changing cores_per_node specific
# to the machine

if machine_name == 'gr':
    machine_long_name = 'grizzly'
    cores_per_node = 36.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    if NODES_REQUIRED < 70:
        qos = "interactive"
    else:
        qos = "standard"
    runcommand = "sbatch -N %d -n %d --qos=%s -J %s -o %s 'performance_test.py' %d %d %s %d %s" % (
        NODES_REQUIRED, max_tasks, qos, job_id, output_name, max_tasks, min_tasks, machine_long_name, cores_per_node, resolution_name)
elif machine_name == 'wf':
    machine_long_name = 'wolf'
    cores_per_node = 16.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    if NODES_REQUIRED < 70:
        qos = "interactive"
    else:
        qos = "standard"
    runcommand = "sbatch -N %d -n %d --qos=%s -J %s -o %s 'performance_test.py' %d %d %s %d %s" % (
        NODES_REQUIRED, max_tasks, qos, job_id, output_name, max_tasks, min_tasks, machine_long_name, cores_per_node, resolution_name)
elif machine_name == 'ba':
    machine_long_name = 'badger'
    cores_per_node = 36.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    if NODES_REQUIRED < 70:
        qos = "interactive"
    else:
        qos = "standard"
    runcommand = "sbatch -N %d -n %d --qos=%s -J %s -o %s 'performance_test.py' %d %d %s %d %s" % (
        NODES_REQUIRED, max_tasks, qos, job_id, output_name, max_tasks, min_tasks, machine_long_name, cores_per_node, resolution_name)
elif machine_name == 'cori-haswell':
    machine_long_name = 'cori-haswell'
    cores_per_node = 32.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    runcommand = "sbatch -N %d -n %d -C haswell --qos=regular -J %s -o %s 'performance_test.py' %d %d %s %d %s" % (
        NODES_REQUIRED, max_tasks, job_id, output_name, max_tasks, min_tasks, machine_long_name, cores_per_node, resolution_name)
elif machine_name == 'cori-knl':
    machine_long_name = 'cori-knl'
    cores_per_node = 68.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    runcommand = "sbatch -N %d -n %d -C knl --qos=regular -J %s -o %s 'performance_test.py' %d %d %s %d %s" % (
        NODES_REQUIRED, max_tasks, job_id, output_name, max_tasks, min_tasks, machine_long_name, cores_per_node, resolution_name)
elif machine_name == 'ed':
    machine_long_name = 'edison'
    cores_per_node = 24.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    runcommand = "sbatch -N %d -n %d --qos=debug -J %s -o %s 'performance_test.py' %d %d %s %d %s" % (
        NODES_REQUIRED, max_tasks, job_id, output_name, max_tasks, min_tasks, machine_long_name, cores_per_node, resolution_name)
elif machine_name[0:5] == 'theta':
    machine_long_name = 'theta'
    cores_per_node = 64.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    runcommand = "qsub -n %d --jobname=%s -O %s 'performance_test.py' %d %d %s %d %s" % (
        NODES_REQUIRED, job_id, output_name, max_tasks, min_tasks, machine_long_name, cores_per_node, resolution_name)
else:
    print "Invalid machine or have not mentioned haswell or knl on Cori"


print "running: ", runcommand
s_args = shlex.split(runcommand)
subprocess.check_call(s_args)
