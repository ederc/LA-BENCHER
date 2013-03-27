/**
 * \file   f4rt.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   February 2013
 * \brief  General source file for f4rt stuff.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "f4rt-config.h"

#include "matrix.h"
#include "dense-mat-mult.h"

#define PACKAGE "F4RT"
#define VERSION "0.0.1"


void print_help(int exval) {
 printf("DESCRIPTION\n");
 printf("       This program is used as a test for parallelization libraries like\n");
 printf("       OpenMP or Intel TBB. It is working with dense floating point entry\n");
 printf("       matrices. Its main features are:\n\n");
 printf("       1. Multiplication of two matrices\n");
 printf("       2. Computing the Gaussian Elimination of a matrix\n\n");
 printf("       Note that when preforming the multiplication he program does not\n");
 printf("       check for overflows. Its purpose is to compare parallelization of mostly\n");
 printf("       independent tasks, namely matrix*matrix kernel operations.\n");
 printf("       The Gaussian Elimination is implemented cache-obliviously as a Gaussian\n");
 printf("       Elimination Paradigm (GEP) as described in \"The Cache-Oblivious Gaussian\n");
 printf("       Elimination Paradigm: Theoretical Framework, Parallelization and Experi-\n");
 printf("       mental Evaluation\" by Chowdhury and Ramachandran.\n\n");
 printf("       The program also provides the possibility of generating random\n");
 printf("       input matrices.\n\n");

 printf("OPTIONS\n");
 printf("       -a        sets CPU affinity of task scheduler; note this only works\n");
 printf("                 with TBB right now\n");
 printf("       -A FILE   set first intput file; if only -A is set but not -B then A*A^T is computed.\n");
 printf("       -b        sets the block size/grain of task scheduler; default = 2\n");
 printf("       -B FILE   set second intput file\n");
 printf("       -C        number of cols of matrix to be generated\n");
 printf("       -d        sets the dimension of the parallel for loop\n");
 printf("                 (values: 1, 2, 3 (3 only for TBB); default = 1)\n");
 printf("       -E        if input file is set, compute the Gaussian Elimination\n");
 printf("       -g        generate a new random float matrix\n");
 printf("       -h        print this help and exit\n");
 printf("       -i        If parallel scheduler is used with option -d1, then the\n");
 printf("                 -i flag triggers collapsing on one loop level deeper.\n");
 printf("                 Without this flag, collapsing is done w.r.t. the outer loop.\n");
 printf("       -m        method to be used: \n");
 printf("                 0 = Sequential\n");
#ifdef __F4RT_HAVE_OPENMP
 printf("                 1 = OpenMP\n");
#endif
#ifdef __F4RT_HAVE_INTEL_TBB
 printf("                 2 = TBB\n");
#endif
#if defined(__F4RT_HAVE_KAAPI)
 printf("                 3 = KAAPI\n");
#endif
#ifdef __F4RT_HAVE_PTHREAD_H
 printf("                 4 = pThread\n");
#endif
 printf("                 Note: By default the sequential implementation is used\n");
 printf("       -M        if input file is set, multiply matrix with its own transpose\n");
 printf("       -N        Not transposing: B is NOT transposed before multiplication,\n");
 printf("                 thus the computation has a way worse cache locality\n");
 printf("       -p        if matrix multiplication took place, print of resulting matrix\n");
 printf("                 (no printing of resulting matrix by default)\n");
 printf("       -R        number of rows of matrix to be generated\n");
 printf("       -s        sets simple task scheduler; note this only works\n");
 printf("                 with TBB right now\n");
 printf("       -t        number of threads to be used (default: all possible ones)\n");
 printf("       -V        print version and exit\n\n");
 printf("       -v        set verbose flag\n");

 exit(exval);
}

void genMatrix(int rows=0, int cols=0) {
  uint32 m, n;
  bool cmp;
  Matrix A;
  if (rows == 0 || cols == 0) {
    std::cout << "Generate new random matrix with entries of type float." << std::endl;
    std::cout << "Number of rows (<2^32): ";
    std::cin >> m;
    std::cout << "Number of cols (<2^32): ";
    std::cin >> n;
    std::cout << "Check if matrix is stored correctly? (1=yes, 0=no)  ";
    std::cin >> cmp;
    A.generateRandomMatrix(m,n,cmp,1);
  } else {
    A.generateRandomMatrix(rows,cols,0,0);
  }
  A.clear();
  std::cout << "Matrix generated." << std::endl;
}

void eliminateMatrix( char* str, int nthrds, int method, int affinity, int blocksize, 
                      int dimension, int impose, int outerloop, int print) {
  Matrix A;

  // read files, stores matrices, etc
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);
  
}

int main(int argc, char *argv[]) {
 int opt;
 char *fileNameA = NULL, *fileNameB = NULL;
 int print = 0, multiply  = 0, nthrds = 0, method = 0, affinity = 0,
     blocksize = 2, dimension = 1, impose = 1, rows = 0, cols = 0,
     generate = 0, outerloop = 1, eliminate = 0;

 /* 
 // no arguments given
 */
 if(argc == 1) {
  //fprintf(stderr, "This program needs arguments....\n\n");
  //print_help(1);
 }

 while((opt = getopt(argc, argv, "hVvgA:B:C:Ept:m:Md:b:aR:Nsi")) != -1) {
  switch(opt) {
    case 'g': 
      generate = 1;
      break;
    case 'h':
      print_help(0);
      break;
    case 'V':
      printf("%s %s\n\n", PACKAGE, VERSION); 
      exit(0);
      break;
    case 'v':
      printf("%s: Verbose option is set `%c'\n", PACKAGE, optopt);
      break;
    case 'A':
      fileNameA  = strdup(optarg);
      //multMatrices(optarg);
      break;
    case 'B':
      fileNameB  = strdup(optarg);
      //multMatrices(optarg);
      break;
    case 'C':
      cols = atoi(strdup(optarg));
      break;
    case 'M':
      multiply  = 1;
      break;
    case 'E':
      eliminate  = 1;
      break;
    case 'p':
      print   = 1;
      break;
    case 'a':
      affinity  = 1;
      break;
    case 's':
      affinity  = 2;
      break;
    case 'N':
      impose  = 0;
      break;
    case 'b':
      blocksize = atoi(strdup(optarg));
      if (blocksize == 0)
        blocksize = 1;
      break;
    case 'd':
      dimension = atoi(strdup(optarg));
      if (dimension == 0)
        dimension = 1;
      if (dimension > 3)
        dimension = 3;
      break;
    case 't':
      nthrds  = atoi(strdup(optarg));
      break;
    case 'm':
      method  = atoi(strdup(optarg));
      break;
    case 'R':
      rows = atoi(strdup(optarg));
      break;
    case 'i':
      outerloop = 0;
      break;
    case ':':
      fprintf(stderr, "%s: Error - Option `%c' needs a value\n\n", PACKAGE, optopt);
      print_help(1);
      break;
    case '?':
      fprintf(stderr, "%s: Error - No such option: `%c'\n\n", PACKAGE, optopt);
      print_help(1);
   }
 }

  /* 
  // print all remaining options
  */
  for(; optind < argc; optind++)
  printf("argument: %s\n", argv[optind]);
   
  if (generate == 1) {
    genMatrix(rows,cols);
  }
  if (multiply && fileNameA) {
    if (fileNameB) {
      multMatrices( fileNameA, fileNameB, nthrds, method, affinity, blocksize,
                    dimension, impose, outerloop,print);  
    } else {
      multEqualMatrices(fileNameA, nthrds, method, affinity, blocksize, 
                        dimension, impose, outerloop, print);  
    }
  }
  if (eliminate && fileNameA) {
      eliminateMatrix(fileNameA, nthrds, method, affinity, blocksize, 
                      dimension, impose, outerloop, print);  
  }
  return 0;
}
