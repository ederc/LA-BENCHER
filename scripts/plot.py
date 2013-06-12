#!/usr/bin/python

import sys
import fnmatch
import os
import glob
import shutil
import argparse
import time
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import pylab as pl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

currentdir = os.getcwd()

parser = argparse.ArgumentParser(description='Visualizes already\
computed benchmarks from F4RT.')
parser.add_argument('-d', '--directory', required=True,
    help='Directory where the benchmark file is included.')
parser.add_argument('-a', '--alg', required=True,
    help='Algorithm used.')
parser.add_argument('-i', '--inc', required=True,
    help='Incremental.')

args = parser.parse_args()

# go to directory
os.chdir(args.directory)

file_name = ''
# find bench file
for files in glob.glob("bench-*"):
  file_name = files

# read lines of the benchmark files
f = open(file_name)
lines = f.readlines()
f.close()
file_name = args.directory
file_name = file_name.replace('/','')
# get
# 1. dimensions of benchmark matrices
# 2. threads for plot, stored in the first line of bench file
dimensions = lines[0].strip().replace(' ','').split(',')
# second line are the thread settings used
plot_threads = lines[1].strip().replace(' ','').split(',')
# get algorithm benchmarked
algorithm_benched = lines[2].strip()
print algorithm_benched
if algorithm_benched == 'Matrix Multiplication':
  algorithm = 1
if algorithm_benched == 'Naive Gaussian Elimination with pivoting':
  algorithm = 2
if algorithm_benched == 'Naive Gaussian Elimination without pivoting':
  algorithm = 3
if algorithm_benched == 'Cache-oblivious Gaussian Elimination without pivoting':
  algorithm = 3
# get threads for plot, stored in the first line of bench file
#plot_threads = f.readline().strip().replace(' ','').split(',')
# for compatibility to the other scripts just store this again
threads = plot_threads
plot_data = list()
start_threads = int(threads[0])
threads = list(map(lambda x: int(x) - 1, threads))

# list of all methods, sequential only if start_threads == 1
if int(args.alg) == 1:
  if start_threads == 1:
    methods = ['Raw sequential','pThread 1D','Open MP collapse(1) outer loop',
    'Open MP collapse(1) inner loop','Open MP collapse(2)',
    'KAAPIC 1D','KAAPIC 2D',
    'Intel TBB 1D auto partitioner','Intel TBB 1D affinity partitioner',
    'Intel TBB 1D simple partitioner','Intel TBB 2D auto partitioner',
    'Intel TBB 2D affinity partitioner','Intel TBB 2D simple partitioner']
  else :
    methods = ['pThread 1D','Open MP collapse(1) outer loop',
    'Open MP collapse(1) inner loop','Open MP collapse(2)',
    'KAAPIC 1D','KAAPIC 2D',
    'Intel TBB 1D auto partitioner','Intel TBB 1D affinity partitioner',
    'Intel TBB 1D simple partitioner','Intel TBB 2D auto partitioner',
    'Intel TBB 2D affinity partitioner','Intel TBB 2D simple partitioner']
if int(args.alg) == 2 or int(args.alg) == 3:
  if start_threads == 1:
    methods = ['Raw sequential','pThread 1D','Open MP collapse(1) outer loop',
    'KAAPIC 1D','Intel TBB 1D auto partitioner','Intel TBB 1D affinity partitioner',
    'Intel TBB 1D simple partitioner']
  else :
    methods = ['pThread 1D','Open MP collapse(1) outer loop',
    'KAAPIC 1D','Intel TBB 1D auto partitioner','Intel TBB 1D affinity partitioner',
    'Intel TBB 1D simple partitioner']
  
if int(args.alg) == 4:
  if start_threads == 1:
    methods = ['Raw sequential','pThread 1D','Open MP parallel sections',
    'KAAPIC Spawn','Intel TBB Invoke']
  else :
    methods = ['pThread 1D','Open MP parallel sections',
    'KAAPIC Spawn','Intel TBB Invoke']

# lists for all methods we have, those are lists of lists:
# E.g. time_series[i] is a list of len(threads) elements of the timings
# of methods[i]. 
time_series = list()
gflops_series = list()
speedup_series = list()

for i in range(0,len(methods)):
  time_series.append(list())
  gflops_series.append(list())
  speedup_series.append(list())


print methods
tmp = -1
tmpold = -1
for l in lines:
  for i in range(0,len(methods)):  
    if l.find(methods[i]) != -1:
      print(methods[i]) 
      tmpold = tmp
      tmp = i
  if l.find('Real time:') != -1:
    if tmpold < tmp:
      time_seq = l.replace('Real time:','').replace('sec','').strip()
      tmpold = tmp
    print(tmp)
    time_series[tmp].append(\
        l.replace('Real time:','').replace('sec','').strip())
    speedup_tmp =float(time_seq) / float(time_series[tmp][len(time_series[tmp])-1])
    print(speedup_tmp)
    speedup_series[tmp].append(str(speedup_tmp))
  if l.find('GFLOPS/sec:') != -1:
    # if the value is inf for infinity due to short computation time, we set
    # the GFLOPS value to be -1
    gflops_series[tmp].append(\
        l.replace('GFLOPS/sec:','').replace('inf','-1').strip())

#plot this data
  #line style, sequential method only if start_threads == 1
  stride = 1
  coloring =\
  [\
  '#0099cc','#33cc00','#ff1b54','#0033cc','#9900cc','#800020',\
  '#ff4011','#ffbf01','#00144f','#ff1450',\
  '#0099cc','#33cc00','#cc0033','#0033cc','#9900cc','#800020',\
  '#ff4011','#ffbf01','#00144f','#ff1450',\
  '#0099cc','#33cc00','#cc0033','#0033cc','#9900cc','#800020',\
  '#ff4011','#ffbf01','#00144f','#ff1450',\
  '#0099cc','#33cc00','#cc0033','#0033cc','#9900cc','#800020',\
  '#ff4011','#ffbf01','#00144f','#ff1450'\
  ]
  styles = [\
  '-','-','-','-','-','-','-','-','-','-',\
  '--','--','--','--','--','--','--','--','--','--',\
  ':',':',':',':',':',':',':',':',':',':',\
  '-','-','-','-','-','-','-','-','-','-'\
  ]
  markers = [\
  'o','None','None','None','None','None','None','None','None','None',\
  'None','None','None','None','None','None','None','None','None','None',\
  'None','None','None','None','None','None','None','None','None','None',\
  'o','o','o','o','o','o','o','o','o','o'\
  ]
  


pl.rc('legend',**{'fontsize':5})
fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('Timings: '+file_name, fontsize=10)
if int(args.alg) == 1:
  pl.title('Mat Mult uint64 Matrix dimensions: '+dimensions[0]+
  ' x '+dimensions[1]+', '+dimensions[1]+' x '+dimensions[2], fontsize=8)
if int(args.alg) == 2 or int(args.alg) == 3:
  if int(args.inc) == -1:
    pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
    ' x '+dimensions[1], fontsize=8)
  else:
    if int(args.inc) == 0:
      pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' with dimensions doubled in each step using '+
      str(max_threads)+' threads', fontsize=8)
    else:
      pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' increasing by '+dimensions[2]+' in each step using '+
      str(max_threads)+' threads', fontsize=8)
if int(args.alg) == 4:
  if int(args.inc) == -1:
    pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
    ' x '+dimensions[1], fontsize=8)
  else:
    if int(args.inc) == 0:
      pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' with dimensions doubled in each step using '+
      str(max_threads)+' threads', fontsize=8)
    else:
      pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' increasing by '+dimensions[2]+' in each step using '+
      str(max_threads)+' threads', fontsize=8)
if int(args.inc) == -1:
  ax.set_xlabel('Number of threads', fontsize=7)
else:
  ax.set_xlabel('Matrix sizes', fontsize=7)
ax.set_ylabel('Real time in seconds', fontsize=8)

pl.grid(b=True, which='major', color='k', linewidth=0.3)
pl.grid(b=True, which='minor', color='k', linewidth=0.1, alpha=0.5)

ax = pl.gca() 

group_labels = plot_threads

#ax.set_xticklabels(group_labels)
threads_tmp = range(0,len(plot_threads))
# get right scale for a4 paper size
scale_tmp = 38 / (len(plot_threads)) 
threads = range(0,38,scale_tmp)
tick_lbs = plot_threads
ax.xaxis.set_ticks(threads)
ax.xaxis.set_ticklabels(tick_lbs)
print(time_series)
print(threads)
p = [None]*len(methods)
for i in range(0,len(methods)):
  p[i], = ax.plot(threads[0:len(time_series[i])], time_series[i], c=coloring[i],
      ls=styles[i], marker=markers[i], markersize='4', markevery=stride, label=i)
# set 0 as min value for y and 1 as min value for x (threads)
#pl.xlim(xmin=1)
pl.ylim(ymin=0)
ax.legend((methods),'upper right', shadow=True, fancybox=True)

# take real time of sequential computation to figure out the 
# granularity of the yaxis
tmp_ticks = ax.yaxis.get_majorticklocs()
granu = tmp_ticks[len(tmp_ticks)-1] / (len(tmp_ticks)-1) / 5
ax.yaxis.set_minor_locator(MultipleLocator(granu))
pl.tick_params(axis='both', which='major', labelsize=6)
pl.tick_params(axis='both', which='minor', labelsize=6)

pl.savefig('timings-plot.pdf',papertype='a4',orientation='landscape')

fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('GFLOPS/sec: '+file_name, fontsize=10)
if int(args.alg) == 1:
  pl.title('Mat Mult uint64 Matrix dimensions: '+dimensions[0]+
  ' x '+dimensions[1]+', '+dimensions[1]+' x '+dimensions[2], fontsize=8)
if int(args.alg) == 2 or int(args.alg) == 3:
  if int(args.inc) == -1:
    pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
    ' x '+dimensions[1], fontsize=8)
  else:
    if int(args.inc) == 0:
      pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' with dimensions doubled in each step using '+
      str(max_threads)+' threads', fontsize=8)
    else:
      pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' increasing by '+dimensions[2]+' in each step using '+
      str(max_threads)+' threads', fontsize=8)
if int(args.alg) == 4:
  if int(args.inc) == -1:
    pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
    ' x '+dimensions[1], fontsize=8)
  else:
    if int(args.inc) == 0:
      pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' with dimensions doubled in each step using '+
      str(max_threads)+' threads', fontsize=8)
    else:
      pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' increasing by '+dimensions[2]+' in each step using '+
      str(max_threads)+' threads', fontsize=8)
if int(args.inc) == -1:
  ax.set_xlabel('Number of threads', fontsize=7)
else:
  ax.set_xlabel('Matrix sizes', fontsize=7)
ax.set_ylabel('GFLOPS per second', fontsize=8)

pl.grid(b=True, which='major', color='k', linewidth=0.3)
pl.grid(b=True, which='minor', color='k', linewidth=0.1, alpha=0.5)

ax = pl.gca() 

#ax.set_xticklabels(group_labels)
threads_tmp = range(0,len(plot_threads))
# get right scale for a4 paper size
scale_tmp = 38 / (len(plot_threads)) 
threads = range(0,38,scale_tmp)
tick_lbs = plot_threads
ax.xaxis.set_ticks(threads)
ax.xaxis.set_ticklabels(tick_lbs)

p = [None]*len(methods)
for i in range(0,len(methods)):
  p[i], = ax.plot(threads[0:len(gflops_series[i])], gflops_series[i], c=coloring[i],
      ls=styles[i], marker=markers[i], markersize='4', label=i)
# set 0 as min value for y and 1 as min value for x (threads)
#pl.xlim(xmin=1)
pl.ylim(ymin=0)
ax.legend((methods),'upper left', shadow=True, fancybox=True)

# take gflops of best computation to figure out the 
# granularity of the yaxis
tmp_ticks = ax.yaxis.get_majorticklocs()
# note that here "abs()" is needed since if computations are too fast we
# set GFLOPS to -1 instead of infinity. Since the MultipleLocator must
# be set to a positive integer value, we have to take care of this case.
granu = abs(tmp_ticks[len(tmp_ticks)-1]) / (len(tmp_ticks)-1) / 5
ax.yaxis.set_minor_locator(MultipleLocator(granu))

pl.tick_params(axis='both', which='major', labelsize=6)
pl.tick_params(axis='both', which='minor', labelsize=6)

pl.savefig('gflops-plot.pdf',papertype='a4',orientation='landscape')

fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('Speedup: '+file_name, fontsize=10)
if int(args.alg) == 1:
  pl.title('Mat Mult uint64 Matrix dimensions: '+dimensions[0]+
  ' x '+dimensions[1]+', '+dimensions[1]+' x '+dimensions[2], fontsize=8)
if int(args.alg) == 2 or int(args.alg) == 3:
  if int(args.inc) == -1:
    pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
    ' x '+dimensions[1], fontsize=8)
  else:
    if int(args.inc) == 0:
      pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' with dimensions doubled in each step using '+
      str(max_threads)+' threads', fontsize=8)
    else:
      pl.title('Naive GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' increasing by '+dimensions[2]+' in each step using '+
      str(max_threads)+' threads', fontsize=8)
if int(args.alg) == 4:
  if int(args.inc) == -1:
    pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
    ' x '+dimensions[1], fontsize=8)
  else:
    if int(args.inc) == 0:
      pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' with dimensions doubled in each step using '+
      str(max_threads)+' threads', fontsize=8)
    else:
      pl.title('Cache-oblivious GEP uint64 Matrix dimensions: '+dimensions[0]+
      ' x '+dimensions[1]+' increasing by '+dimensions[2]+' in each step using '+
      str(max_threads)+' threads', fontsize=8)
if int(args.inc) == -1:
  ax.set_xlabel('Number of threads', fontsize=7)
else:
  ax.set_xlabel('Matrix sizes', fontsize=7)
ax.set_ylabel('Speedup', fontsize=8)

pl.grid(b=True, which='major', color='k', linewidth=0.3)
pl.grid(b=True, which='minor', color='k', linewidth=0.1, alpha=0.5)

ax = pl.gca() 

#ax.set_xticklabels(group_labels)
threads_tmp = range(0,len(plot_threads))
# get right scale for a4 paper size
scale_tmp = 38 / (len(plot_threads)) 
threads = range(0,38,scale_tmp)
tick_lbs = plot_threads
ax.xaxis.set_ticks(threads)
ax.xaxis.set_ticklabels(tick_lbs)

p = [None]*len(methods)
for i in range(0,len(methods)):
  p[i], = ax.plot(threads[0:len(speedup_series[i])], speedup_series[i], c=coloring[i],
      ls=styles[i], marker=markers[i], markersize='4', label=i)
# set 0 as min value for y and 1 as min value for x (threads)
#pl.xlim(xmin=1)
pl.ylim(ymin=0)
ax.legend((methods),'upper left', shadow=True, fancybox=True)

# take gflops of best computation to figure out the 
# granularity of the yaxis
tmp_ticks = ax.yaxis.get_majorticklocs()
# note that here "abs()" is needed since if computations are too fast we
# set GFLOPS to -1 instead of infinity. Since the MultipleLocator must
# be set to a positive integer value, we have to take care of this case.
granu = abs(tmp_ticks[len(tmp_ticks)-1]) / (len(tmp_ticks)-1) / 5
ax.yaxis.set_minor_locator(MultipleLocator(granu))

pl.tick_params(axis='both', which='major', labelsize=6)
pl.tick_params(axis='both', which='minor', labelsize=6)

pl.savefig('speedup-plot.pdf',papertype='a4',orientation='landscape')
