#!/usr/bin/python

import sys
import fnmatch
import os
import shutil
import argparse
import time
import math
import pylab as pl
from matplotlib.ticker import multiplelocator, formatstrformatter

# generates a 10 digit hash value from the current date and time
# to append to the already existing ".singularrc" filename in order
# to restore it later on
def hash():
  return '{0:010x}'.format(int(time.time() * 256))[:10]

# gets number of lines in a given text file
def linesinfile(file_name):
    with open(file_name) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
  
currentdir = os.getcwd()

parser = argparse.argumentparser(description='generates two random matrices\
a and b with float entries. these matrices are then multiplied, thus it\
is enough to predefine #rows of a, #cols of a and #cols of b. the\
multiplication is performed in various combinations of the f4rt dense\
matrix multiplication implementation. afterwards the stored results\
are visualized.')
parser.add_argument('-l', '--rowsa', required=true,
    help='number of rows of matrix a')
parser.add_argument('-m', '--colsa', required=true,
    help='number of cols of matrix a')
parser.add_argument('-n', '--colsb', required=true,
    help='number of cols of matrix b')
parser.add_argument('-t', '--threads', required=true,
    help='maximal number of threads to be used')
parser.add_argument('-b', '--base', default=2,
    help='base of number of threads, e.g. -b 2 -t 16 would lead to computations\
in 1,2,4,8,16 threads. Default is 2.')
parser.add_argument('-hg', '--homog', action="store_true", default=false,
    help='set this flag if the ideal should be homogenized.\
note that this flag increases the number of variables by 1.')

args = parser.parse_args()

# range of threads
threads = list()
exp = 0
max_threads = int(args.threads)
base = int(args.base)
if base == 1:
  while (base) <= max_threads:
    threads.append(base)
    base += 1
else:
  while (base**exp) <= max_threads:
    threads.append(base**exp)
    exp += 1

# list of all methods
methods = ['raw sequential','open mp collapse(1)','open mp collapse(2)',
'intel tbb 1d auto partitioner','intel tbb 1d affinity partitioner',
'intel tbb 1d simple partitioner','intel tbb 2d auto partitioner',
'intel tbb 2d affinity partitioner','intel tbb 2d simple partitioner']

# lists for all methods we have, those are lists of lists:
# e.g. time_series[i] is a list of len(threads) elements of the timings
# of methods[i]. 
time_series = list()
gflops_series = list()

for i in range(0,len(methods)):
  time_series.append(list())
  gflops_series.append(list())

# generate hash value if needed
hash_value = os.urandom(16).encode('hex')

folder_name = "test-"+str(hash_value)

if not os.path.exists(folder_name):
  os.makedirs(folder_name)

os.chdir(os.getcwd()+"/"+folder_name)

#generate random matrices without timestamp
os.system('../../src/dense-mult -g -r '+args.rowsa+' -c '+args.colsa)
os.system('../../src/dense-mult -g -r '+args.colsa+' -c '+args.colsb)

bench_file = "bench-"+str(hash_value)
f = open(bench_file,"w")

strstr = '../../src/dense-mult -c \
-a random-float-mat-'+args.rowsa+'-'+args.colsa+'.mat \
-b random-float-mat-'+args.colsa+'-'+args.colsb+'.mat'

thrds_str = str(threads)
thrds_str = thrds_str.replace('[','')
thrds_str = thrds_str.replace(']','')
thrds_str = thrds_str
f.write(thrds_str+'\n')
f.close()

# sequential computation
print(strstr+' -i -m0 >> '+bench_file+'...')
os.system(strstr+' -i -m0 >> '+bench_file)
print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# OpenMP computations 1D
for i in threads:
  print(strstr+' -i -m1 -d 1 -t '+str(i)+' >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m1 -d 1 -t '+str(i)+' >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# OpenMP computations 2D
for i in threads:
  print(strstr+' -i -m1 -d 2 -t '+str(i)+' >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m1 -d 2 -t '+str(i)+' >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D auto
for i in threads:
  print(strstr+' -i -m2 -t '+str(i)+' >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m2 -t '+str(i)+' >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D affinity
for i in threads:
  print(strstr+' -i -m2 -t '+str(i)+' -a >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m2 -t '+str(i)+' -a >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D simple
for i in threads:
  print(strstr+' -i -m2 -t '+str(i)+' -s >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m2 -t '+str(i)+' -s >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D auto
for i in threads:
  print(strstr+' -i -m2 -t '+str(i)+' -d 2 >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m2 -t '+str(i)+' -d 2 >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D affinity
for i in threads:
  print(strstr+' -i -m2 -t '+str(i)+' -d 2 -a >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m2 -t '+str(i)+' -d 2 -a >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D simple
for i in threads:
  print(strstr+' -i -m2 -t '+str(i)+' -d 2 -s >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m2 -t '+str(i)+' -d 2 -s >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# read lines of the benchmark files
f = open('bench-'+str(hash_value))

# get threads for plot, stored in the first line of bench file
plot_threads = f.readline().strip().replace(' ','').split(',')
# for compatibility to the other scripts just store this again
threads = plot_threads
threads = list(map(lambda x: int(x) - 1, threads))
lines = f.readlines()
f.close()


for l in lines:
  for i in range(0,len(methods)):  
    if l.find(methods[i]) != -1:
      tmp = i
  if l.find('real time:') != -1:
    time_series[tmp].append(\
        l.replace('real time:','').replace('sec','').strip())
  if l.find('gflops/sec:') != -1:
    # if the value is inf for infinity due to short computation time, we set
    # the gflops value to be -1
    gflops_series[tmp].append(\
        l.replace('gflops/sec:','').replace('inf','-1').strip())

#plot this data

#line style
coloring = ['k^','b-','b--','g-','g--','g:','r-','r--','r:']

pl.rc('legend',**{'fontsize':5})
fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('timings: '+file_name, fontsize=10)
ax.set_xlabel('number of threads', fontsize=8)
ax.set_ylabel('real time in seconds', fontsize=8)

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

p = [none]*len(methods)
for i in range(0,len(methods)):
  p[i], = ax.plot(threads[0:len(time_series[i])], time_series[i], coloring[i], label=i)
# set 0 as min value for y and 1 as min value for x (threads)
#pl.xlim(xmin=1)
pl.ylim(ymin=0)
ax.legend((methods),'upper right', shadow=true, fancybox=true)

# take real time of sequential computation to figure out the 
# granularity of the yaxis
tmp_ticks = ax.yaxis.get_majorticklocs()
granu = tmp_ticks[len(tmp_ticks)-1] / (len(tmp_ticks)-1) / 5
ax.yaxis.set_minor_locator(multiplelocator(granu))
pl.tick_params(axis='both', which='major', labelsize=6)
pl.tick_params(axis='both', which='minor', labelsize=6)

pl.savefig('timings-plot.pdf',papertype='a4',orientation='landscape')

fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('gflops/sec: '+file_name, fontsize=10)
ax.set_xlabel('number of threads', fontsize=8)
ax.set_ylabel('gflops per second', fontsize=8)

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

p = [none]*len(methods)
for i in range(0,len(methods)):
  p[i], = ax.plot(threads[0:len(gflops_series[i])], gflops_series[i],coloring[i], label=i)
# set 0 as min value for y and 1 as min value for x (threads)
#pl.xlim(xmin=1)
pl.ylim(ymin=0)
ax.legend((methods),'upper left', shadow=true, fancybox=true)

# take gflops of best computation to figure out the 
# granularity of the yaxis
tmp_ticks = ax.yaxis.get_majorticklocs()
# note that here "abs()" is needed since if computations are too fast we
# set gflops to -1 instead of infinity. since the multiplelocator must
# be set to a positive integer value, we have to take care of this case.
granu = abs(tmp_ticks[len(tmp_ticks)-1]) / (len(tmp_ticks)-1) / 5
ax.yaxis.set_minor_locator(multiplelocator(granu))

pl.tick_params(axis='both', which='major', labelsize=6)
pl.tick_params(axis='both', which='minor', labelsize=6)

pl.savefig('gflops-plot.pdf',papertype='a4',orientation='landscape')
