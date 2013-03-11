#!/usr/bin/python

import sys
import fnmatch
import os
import shutil
import argparse
import time

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

# generate hash value if needed
hash_value = hash()

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
print(strstr+' -m 2 >> bench-'+str(hash_value)+'...')
os.system(strstr+' -m 2 >> bench-'+str(hash_value))
print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# OpenMP computations
for i in range(0,6):
  print(strstr+' -m 1 -t '+str(2**i)+' >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -m 1 -t '+str(2**i)+' >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D auto
for i in range(0,6):
  print(strstr+' -m 0 -t '+str(2**i)+' >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -m 0 -t '+str(2**i)+' >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D affinity
for i in range(0,6):
  print(strstr+' -m 0 -t '+str(2**i)+' -a >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -m 0 -t '+str(2**i)+' -a >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 1D simple
for i in range(0,6):
  print(strstr+' -m 0 -t '+str(2**i)+' -s >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -m 0 -t '+str(2**i)+' -s >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D auto
for i in range(0,6):
  print(strstr+' -m 0 -t '+str(2**i)+' -d 2 >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -m 0 -t '+str(2**i)+' -d 2 >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D affinity
for i in range(0,6):
  print(strstr+' -m 0 -t '+str(2**i)+' -d 2 -a >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -m 0 -t '+str(2**i)+' -d 2 -a >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

# TBB computations 2D simple
for i in range(0,6):
  print(strstr+' -m 0 -t '+str(2**i)+' -d 2 -s >> bench-'+str(hash_value)+'...')
  os.system(strstr+' -m 0 -t '+str(2**i)+' -d 2 -s >> bench-'+str(hash_value))
  print 'Done at '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())
