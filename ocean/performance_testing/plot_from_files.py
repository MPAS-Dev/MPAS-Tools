#!/usr/bin/env python

"""
Name: plot_from_files.py
Author: Divya Jaganathan
Date: 26 July, 2018

Plots a single plot of different performance curves from different performance_data text files in a folder

"""

import glob
import matplotlib
matplotlib.use('agg')
import numpy as np
import matplotlib.pyplot as plt
import subprocess

path = "/lustre/scratch2/turquoise/divjag2005/case_runs/performance_results_4096/*.txt"

files = glob.glob(path)
files = [files[4], files[0], files[1], files[3], files[2]]
num_files = len(files)
print(num_files)

no_res_in_a_file = 9
file_counter = 0

array_x = np.zeros(shape=(num_files, no_res_in_a_file))
array_y = np.zeros(shape=(num_files, no_res_in_a_file))

colors = ["g", "k", "m", "r", "b"]
labels = [
    "uniform 60km",
    "variable 60to30km",
    "uniform 30km",
    "variable 60to15km",
    "uniform 15km"]


for file in files:

    f = open(file, 'r')
    ob = f.read().split('\n')
    num_lines = len(ob) - 1
    line_counter = 6
    i = 0
    rank_column = 0
    SYPD_column = 7

    while line_counter < num_lines:
        array_x[file_counter][i] = ob[line_counter].split('\t')[rank_column]
        array_y[file_counter][i] = ob[line_counter].split('\t')[SYPD_column]
        line_counter = line_counter + 1
        i = i + 1

    font = {'weight': 'bold',
            'size': '14'}

    matplotlib.rc('font', **font)
    plt.loglog(array_x[file_counter][0:i -
                                     1], array_y[file_counter][0:i -
                                                               1], '-o', color=colors[file_counter], label="%s" %
               labels[file_counter])
    perfect = (array_y[file_counter][i - 1] /
               array_x[file_counter][i - 1]) * array_x[file_counter][0:i - 1]
    plt.loglog(
        array_x[file_counter][0:i - 1],
        perfect,
        '--',
        color=colors[file_counter])
    file_counter = file_counter + 1
    f.close()

plt.xlabel('Number of MPI ranks', fontsize=14, weight='bold')
plt.ylabel('SYPD', fontsize=14, weight='bold')
plt.title(' 36 Core Broadwell (No HT)', fontsize=14, weight='bold')
plt.xlim((10, 10000))
plt.ylim((0.05, 4000))
plt.tight_layout()
plt.grid()
plt.legend(title='resolution', loc='upper left')
plt.savefig('result.png')
