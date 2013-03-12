#!/usr/bin/python

import sys
import fnmatch
import os
import shutil
import argparse
import time
import pylab as pl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# generates a 10 digit hash value from the current date and time
# to append to the already existing ".singularrc" filename in order
# to restore it later on
def hash():
  return '{0:010x}'.format(int(time.time() * 256))[:10]

# gets number of lines in a given text file
def linesInFile(file_name):
    with open(file_name) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
  
currentdir = os.getcwd()

parser = argparse.ArgumentParser(description='Generates two random matrices\
A and B with float entries. These matrices are then multiplied, thus it\
is enough to predefine #rows of A, #cols of A and #cols of B. The\
multiplication is performed in various combinations of the F4RT dense\
matrix multiplication implementation. Afterwards the stored results\
are visualized.')
parser.add_argument('-l', '--rowsa', required=True,
    help='Number of rows of matrix A')
parser.add_argument('-m', '--colsa', required=True,
    help='Number of cols of matrix A')
parser.add_argument('-n', '--colsb', required=True,
    help='Number of cols of matrix B')
parser.add_argument('-hg', '--homog', action="store_true", default=False,
    help='Set this flag if the ideal should be homogenized.\
    Note that this flag increases the number of variables by 1.')

args = parser.parse_args()

# range of threads
threads = range(0,6)

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

# generate hash value if needed
hash_value = os.urandom(16).encode('hex')

folder_name = "test-"+str(hash_value)

if not os.path.exists(folder_name):
  os.makedirs(folder_name)

os.chdir(os.getcwd()+"/"+folder_name)

#generate random matrices without timestamp
os.system('../../src/dense-mult -g -R '+args.rowsa+' -C '+args.colsa)
os.system('../../src/dense-mult -g -R '+args.colsa+' -C '+args.colsb)

#f = open("bench-"+str(hash_value),"a")

strstr = '../../src/dense-mult -c \
-A random-float-mat-'+args.rowsa+'-'+args.colsa+'.mat \
-B random-float-mat-'+args.colsa+'-'+args.colsb+'.mat'

# sequential computation
print(strstr+' -i -m 2 >> bench-'+str(hash_value)+'...')
os.system(strstr+' -i -m 2 >> bench-'+str(hash_value))
print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# OpenMP computations
for i in threads:
  print(strstr+' -i -m 1 -t '+str(2**i)+' >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m 1 -t '+str(2**i)+' >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D auto
for i in threads:
  print(strstr+' -i -m 0 -t '+str(2**i)+' >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m 0 -t '+str(2**i)+' >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D affinity
for i in threads:
  print(strstr+' -i -m 0 -t '+str(2**i)+' -a >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m 0 -t '+str(2**i)+' -a >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D simple
for i in threads:
  print(strstr+' -i -m 0 -t '+str(2**i)+' -s >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m 0 -t '+str(2**i)+' -s >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D auto
for i in threads:
  print(strstr+' -i -m 0 -t '+str(2**i)+' -d 2 >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m 0 -t '+str(2**i)+' -d 2 >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D affinity
for i in threads:
  print(strstr+' -i -m 0 -t '+str(2**i)+' -d 2 -a >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m 0 -t '+str(2**i)+' -d 2 -a >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D simple
for i in threads:
  print(strstr+' -i -m 0 -t '+str(2**i)+' -d 2 -s >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -i -m 0 -t '+str(2**i)+' -d 2 -s >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# read lines of the benchmark files
f = open('bench-'+str(hash_value))
lines = f.readlines()
f.close()

for l in lines:
  for i in range(0,len(methods)):  
    if l.find(methods[i]) != -1:
      print methods[i]
      print i
      tmp = i
  if l.find('Real time:') != -1:
    print 'real time'
    time_series[tmp].append(\
        l.replace('Real time:','').replace('sec','').strip())
  if l.find('GFLOPS/sec:') != -1:
    print 'gflops'
    # if the value is inf for infinity due to short computation time, we set
    # the GFLOPS value to be -1
    gflops_series[tmp].append(\
        l.replace('GFLOPS/sec:','').replace('inf','-1').strip())

print time_series
print gflops_series

#plot this data

pl.rc('legend',**{'fontsize':6})
fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('Timing results', fontsize=14, fontweight='bold')
ax.set_xlabel('Number of threads')
ax.set_ylabel('Real time in seconds')

#pl.grid(True)

ax = pl.gca() 

ax.xaxis.set_minor_locator(MultipleLocator(0.25))

group_labels = ['1','2','4','8','16','32']

ax.set_xticklabels(group_labels)

p = [None]*len(methods)
for i in range(0,len(methods)):
  p[i], = ax.plot(threads[0:len(time_series[i])], time_series[i],coloring[i], label=i)
ax.legend((methods),'upper right', shadow=True, fancybox=True)

pl.savefig('timings-plot.pdf')

fig = pl.figure()
ax = fig.add_subplot(111)
fig.suptitle('GFLOPS results', fontsize=14, fontweight='bold')
ax.set_xlabel('Number of threads')
ax.set_ylabel('GFLOPS per seconds')

#pl.grid(True)

ax = pl.gca() 

ax.xaxis.set_minor_locator(MultipleLocator(0.25))

group_labels = ['1','2','4','8','16','32']

ax.set_xticklabels(group_labels)

p = [None]*len(methods)
for i in range(0,len(methods)):
  p[i], = ax.plot(threads[0:len(gflops_series[i])], gflops_series[i],coloring[i], label=i)
ax.legend((methods),'lower right', shadow=True, fancybox=True)

pl.savefig('gflops-plot.pdf')
