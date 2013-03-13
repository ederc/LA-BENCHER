#!/usr/bin/python

import sys
import fnmatch
import os
import glob
import shutil
import argparse
import time
import pylab as pl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

currentdir = os.getcwd()

parser = argparse.ArgumentParser(description='Visualizes already\
computed benchmarks from F4RT.')
parser.add_argument('-d', '--directory', required=True,
    help='Directory where the benchmark file is included.')

args = parser.parse_args()

# list of all methods
methods = ['Raw sequential','Open MP collapse(1)','Open MP collapse(2)',
'Intel TBB 1D auto partitioner','Intel TBB 1D affinity partitioner',
'Intel TBB 1D simple partitioner','Intel TBB 2D auto partitioner',
'Intel TBB 2D affinity partitioner','Intel TBB 2D simple partitioner']

coloring = ['k-','b--','b-.','g--','g-.','g:','r--','r-.','r:']
# lists for all methods we have, those are lists of lists:
# E.g. time_series[i] is a list of len(threads) elements of the timings
# of methods[i]. 
time_series = list()
gflops_series = list()

for i in range(0,len(methods)):
  time_series.append(list())
  gflops_series.append(list())

# go to directory
os.chdir(args.directory)

file_name = ''
# find bench file
for files in glob.glob("bench-*"):
  file_name = files

# read lines of the benchmark files
f = open(file_name)

# get threads for plot, stored in the first line of bench file
plot_threads = f.readline().rstrip().split(',')
# for compatibility to the other scripts just store this again
threads = plot_threads
lines = f.readlines()
f.close()


for l in lines:
  for i in range(0,len(methods)):  
    if l.find(methods[i]) != -1:
      tmp = i
  if l.find('Real time:') != -1:
    time_series[tmp].append(\
        l.replace('Real time:','').replace('sec','').strip())
  if l.find('GFLOPS/sec:') != -1:
    # if the value is inf for infinity due to short computation time, we set
    # the GFLOPS value to be -1
    gflops_series[tmp].append(\
        l.replace('GFLOPS/sec:','').replace('inf','-1').strip())

#plot this data

pl.rc('legend',**{'fontsize':6})
fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('Timing results', fontsize=14, fontweight='bold')
ax.set_xlabel('Number of threads')
ax.set_ylabel('Real time in seconds')

#pl.grid(True)

ax = pl.gca() 

group_labels = plot_threads

ax.set_xticklabels(group_labels)

p = [None]*len(methods)
for i in range(0,len(methods)):
  print i
  print time_series[i]
  #p[i], = ax.plot(threads[0:len(time_series[i])], time_series[i], coloring[i], label=i)
ax.legend((methods),'upper right', shadow=True, fancybox=True)

ax.xaxis.set_minor_locator(MultipleLocator(0.20))

# take real time of sequential computation to figure out the 
# granularity of the yaxis
tmp_ticks = ax.yaxis.get_majorticklocs()
granu = tmp_ticks[len(tmp_ticks)-1] / (len(tmp_ticks)-1) / 5
ax.yaxis.set_minor_locator(MultipleLocator(granu))

pl.savefig('timings-plot.pdf')

fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('GFLOPS results', fontsize=14, fontweight='bold')
ax.set_xlabel('Number of threads')
ax.set_ylabel('GFLOPS per seconds')

ax = pl.gca() 

ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

group_labels = plot_threads

ax.set_xticklabels(group_labels)

p = [None]*len(methods)
for i in range(0,len(methods)):
  p[i], = ax.plot(threads[0:len(gflops_series[i])], gflops_series[i],coloring[i], label=i)
ax.legend((methods),'upper left', shadow=True, fancybox=True)

# take gflops of best computation to figure out the 
# granularity of the yaxis
tmp_ticks = ax.yaxis.get_majorticklocs()
# note that here "abs()" is needed since if computations are too fast we
# set GFLOPS to -1 instead of infinity. Since the MultipleLocator must
# be set to a positive integer value, we have to take care of this case.
granu = abs(tmp_ticks[len(tmp_ticks)-1]) / (len(tmp_ticks)-1) / 5
ax.yaxis.set_minor_locator(MultipleLocator(granu))

pl.savefig('gflops-plot.pdf')
