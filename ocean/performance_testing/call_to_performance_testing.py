#!/usr/bin/env python

"""
Name: call_to_performance_testing.py
Author: Divya Jaganathan
Date: July 6, 2018

command format: python call_to_performance_testing.py -c <machine> -M <Maximum Tasks> -m <Minimum Tasks (Optional,default=2)> -p <plot_only=true/false>

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
    "-c",
    "--cpu_type",
    dest="cpu_type",
    help="If cori, enter cori-haswell/cori-knl",
    default=os.uname()[1][0:2])
parser.add_argument(
    "-M",
    "--max_tasks",
    dest="max_tasks",
    help="Maximum number of tasks",
    required=True)
parser.add_argument(
    "-m",
    "--min_tasks",
    dest="min_tasks",
    help="Minimum number of tasks",
    default=2)
parser.add_argument(
    "-r",
    "--resolution",
    dest="resolution",
    help="Resolution ",
    default="QU")
args = parser.parse_args()

cpu_type = args.cpu_type
max_tasks = int(args.max_tasks)
min_tasks = int(args.min_tasks)
res = args.resolution

job_id = res + "_perf_" + str(max_tasks)
output_name = "slurm_" + job_id + ".out"

if cpu_type == 'gr':
    cores_per_node = 36.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    if NODES_REQUIRED < 70:
        qos = "interactive"
    else:
        qos = "standard"
    runcommand = "sbatch -N %d -n %d --qos=%s -J %s -o %s 'performance_testing.py' %d %d %s %d" % (
        NODES_REQUIRED, max_tasks, qos, job_id, output_name, max_tasks, min_tasks, cpu_type, cores_per_node)
elif cpu_type == 'cori-haswell':
    cores_per_node = 32.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    runcommand = "sbatch -N %d -n %d -C haswell --qos=regular -J %s -o %s 'performance_testing.py' %d %d %s %d" % (
        NODES_REQUIRED, max_tasks, job_id, output_name, max_tasks, min_tasks, cpu_type, cores_per_node)
elif cpu_type == 'cori-knl':
    cores_per_node = 68.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    runcommand = "sbatch -N %d -n %d -C knl --qos=regular -J %s -o %s 'performance_testing.py' %d %d %s %d" % (
        NODES_REQUIRED, max_tasks, job_id, output_name, max_tasks, min_tasks, cpu_type, cores_per_node)
elif cpu_type == 'ed':
    cores_per_node = 24.0
    NODES_REQUIRED = int(np.ceil(max_tasks / cores_per_node))
    runcommand = "sbatch -N %d -n %d --qos=debug -J %s -o %s 'performance_testing.py' %d %d %s %d" % (
        NODES_REQUIRED, max_tasks, job_id, output_name, max_tasks, min_tasks, cpu_type, cores_per_node)
else:
    print "Invalid Machine or have not mentioned haswell or knl on Cori"


s_args = shlex.split(runcommand)
print "running", ''.join(s_args)

subprocess.check_call(s_args)
