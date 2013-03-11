#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "matrix-tbb.h"
#include "matrix-omp.h"
#include "matrix-seq.h"

#define PACKAGE "dense-mult"
#define VERSION "0.0.1"


void print_help(int exval) {
 printf("DESCRIPTION\n");
 printf("       This program is used as a test for parallelization libraries like\n");
 printf("       OpenMP or Intel TBB. Its main feature is to multiply 2 dense matrices\n");
 printf("       with uint32 entries. If the value of the entries are < 2^16 the\n");
 printf("       result is even correct. Otherwise the program does not check for\n");
 printf("       overflows. Its purpose is to compare parallelization of mostly\n");
 printf("       independent tasks, namely matrix*matrix kernel operations.\n\n");
 printf("       The program also provides the possibility of generating random\n");
 printf("       input matrices.\n\n");

 printf("OPTIONS\n");
 printf("       -A FILE   set first intput file; if only -A is set but not -B then A*A^T is computed.\n");
 printf("       -B FILE   set second intput file\n");
 printf("       -a        sets CPU affinity of task scheduler; note this only works\n");
 printf("                 with TBB right now\n");
 printf("       -b        sets the block size/grain of task scheduler; default = 2\n");
 printf("       -c        if input file is set, multiply matrix with its own transpose\n");
 printf("       -d        sets the dimension of the parallel for loop; note this only works\n");
 printf("                 with TBB right now (possible values: 1, 2; default = 1)\n");
 printf("       -g        generate a new random uint16 matrix\n");
 printf("       -h        print this help and exit\n");
 printf("       -i        impose, i.e. cheat: Transpose B before multiplication and use\n");
 printf("                 better cache locality\n");
 printf("       -m        method to be used: \n");
 printf("                 0 = TBB,\n");
 printf("                 1 = OpenMP\n");
 printf("                 2 = Sequential\n");
 printf("                 Note: By default TBB is used\n");
 printf("       -p        if matrix multiplication took place, print of resulting matrix\n");
 printf("                 (no printing of resulting matrix by default)\n");
 printf("       -s        sets simple task scheduler; note this only works\n");
 printf("                 with TBB right now\n");
 printf("       -t        number of threads to be used (default: all possible ones)\n");
 printf("       -V        print version and exit\n\n");
 printf("       -v        set verbose flag\n");

 exit(exval);
}

void genMatrix() {
  uint32 m, n;
  bool cmp;
  std::cout << "Generate new random matrix with entries of type float." << std::endl;
  std::cout << "Number of rows (<2^32): ";
  std::cin >> m;
  std::cout << "Number of cols (<2^32): ";
  std::cin >> n;
  Matrix A;
  std::cout << "Check if matrix is stored correctly? (1=yes, 0=no)  ";
  std::cin >> cmp;
  A.generateRandomMatrix(m,n,cmp);
  A.clear();
  std::cout << "Matrix generated." << std::endl;
}

void prepareMult(Matrix& A, Matrix& B, char* str) {
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);

  // let B be just a copy of A
  // we will then multiply A*B^T
  B.transpose(A);
}

void multiply(Matrix& C, const Matrix& A, const Matrix& B, const int nthrds, const int blocksize,
              const int method, const int dimension, const int affinity, int impose) {
  // C = A*B^T
  if (method == 0) { // TBB
    if (dimension == 1) {
      if (affinity == 1) {
        multTBBAffine(C, A, B, nthrds, blocksize, impose);
      } else {
        if (affinity == 2) {
          multTBBSimple(C, A, B, nthrds, blocksize, impose);
        } else {
          multTBBAuto(C, A, B, nthrds, blocksize, impose);
        }
      }
    } else {
      if (affinity == 1) {
        multTBBAffine2d(C, A, B, nthrds, blocksize, impose);
      } else {
        if (affinity == 2) {
          multTBBSimple2d(C, A, B, nthrds, blocksize, impose);
        } else {
          multTBBAuto2d(C, A, B, nthrds, blocksize, impose);
        }
      }
    }
  }
  if (method == 1) // OpenMP
    multOMP(C, A, B, nthrds, blocksize, impose);
  if (method == 2) // plain sequential w/o scheduler overhead
    multSEQ(C, A, B, blocksize, impose);
}

void multMatrices(char* str1, char* str2, int nthrds, int method, int affinity, int blocksize, 
                  int dimension, int impose, int print) {
  Matrix A, B;

  // read files, stores matrices, etc
  FILE* file1  = fopen(str1,"rb");
  A.read(file1);
  FILE* file2  = fopen(str2,"rb");
  B.read(file2);


  if (impose == 1) {
    B.transpose();
    // check dimensions
    if (A.nCols() != B.nCols()) {
      fprintf(stderr, "Dimensions of A and B are not correct!\nProgram exiting now!\n");
      exit(EXIT_FAILURE);
    }
  } else {
    // check dimensions
    if (A.nCols() != B.nRows()) {
      fprintf(stderr, "Dimensions of A and B are not correct!\nProgram exiting now!\n");
      exit(EXIT_FAILURE);
    }
  }
  
  //B.print();
  Matrix C(A.nRows(), B.nCols());

  multiply(C, A, B, nthrds, blocksize, method, dimension, affinity, impose);

  if (print)
    C.print();
  // clear memory
  A.clear();
  B.clear();
  C.clear();
}

void multEqualMatrices( char* str, int nthrds, int method, int affinity, int blocksize, 
                        int dimension, int impose, int print) {
  Matrix A, B;

  // read files, stores matrices, etc
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);

  // let B be just a copy of A
  // we will then multiply A*B^T
  B.transpose(A);
  
  Matrix C(A.nRows(), B.nRows());

  multiply(C, A, B, nthrds, blocksize, method, dimension, affinity,impose);

  if (print)
    C.print();
  // clear memory
  A.clear();
  B.clear();
  C.clear();
}

int main(int argc, char *argv[]) {
 int opt;
 char *fileNameA = NULL, *fileNameB = NULL;
 int print = 0, multiply  = 0, nthrds = 0, method = 0, affinity = 0,
     blocksize = 2, dimension = 1, impose = 0;

 /* 
 // no arguments given
 */
 if(argc == 1) {
  //fprintf(stderr, "This program needs arguments....\n\n");
  //print_help(1);
 }

 while((opt = getopt(argc, argv, "hVvgA:B:pt:m:cd:b:aiso:")) != -1) {
  switch(opt) {
    case 'g': 
      genMatrix();
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
    case 'c':
      multiply  = 1;
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
    case 'i':
      impose  = 1;
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
      if (dimension > 2)
        dimension = 2;
      break;
    case 't':
      nthrds  = atoi(strdup(optarg));
      break;
    case 'm':
      method  = atoi(strdup(optarg));
      break;
    case 'o':
      printf("Output: %s\n", optarg);
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

  if (multiply && fileNameA) {
    if (fileNameB) {
      multMatrices(fileNameA, fileNameB, nthrds, method, affinity, blocksize, dimension, impose, print);  
    } else {
      multEqualMatrices(fileNameA, nthrds, method, affinity, blocksize, dimension, impose, print);  
    }
  }
 return 0;
}
