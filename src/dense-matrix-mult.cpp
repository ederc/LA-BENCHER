#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "matrix-tbb.h"
#include "matrix-omp.h"

#define PACKAGE "dense-mult"
#define VERSION "0.0.1"


void print_help(int exval) {
 printf("%s,%s show working getopt example\n", PACKAGE, VERSION); 
 printf("%s [-h] [-V] [-f FILE] [-o FILE]\n\n", PACKAGE);

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
 printf("       -h        print this help and exit\n");
 printf("       -V        print version and exit\n\n");

 printf("       -v        set verbose flag\n");
 printf("       -g        generate a new random uint16 matrix\n");
 printf("       -f FILE   set intput file\n");
 printf("       -c        if input file is set, multiply matrix with its own transpose\n");
 printf("       -p        if matrix multiplication took place, print of resulting matrix\n");
 printf("                 (no printing of resulting matrix by default)\n");
 printf("       -t        number of threads to be used (default: all possible ones)\n");
 printf("       -m        method to be used: 0=TBB, 1=OpenMP (default: TBB)\n");

 exit(exval);
}

void genMatrix() {
  uint32 m, n;
  bool cmp;
  std::cout << "Generate new random matrix with entries of type uint16." << std::endl;
  std::cout << "Number of rows (<4294967296): ";
  std::cin >> m;
  std::cout << "Number of cols (<4294967296): ";
  std::cin >> n;
  Matrix A;
  std::cout << "Check if matrix is stored correctly? (1=yes, 0=no)  ";
  std::cin >> cmp;
  A.generateRandomMatrix(m,n,cmp);
  std::cout << "Matrix generated." << std::endl;
}

void prepareMult(Matrix& A, Matrix& B, char* str) {
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);

  // let B be just a copy of A
  // we will then multiply A*B^T
  B.copy(A);
}

void multMatrices(char* str, int nthrds, int method, int print) {
  Matrix A, B;

  // read files, stores matrices, etc
  prepareMult(A, B, str);
  
  Matrix C(A.nRows(), B.nRows());

  // C = A*B^T
  if (!method)
    multTBB(C, A, B, nthrds);
  else
    multOMP(C, A, B, nthrds);
  if (print)
    C.print();
  // clear memory
  A.clear();
  B.clear();
  C.clear();
}

int main(int argc, char *argv[]) {
 int opt;
 char* fileName;
 int print = 0, multiply  = 0, nthrds = 0, method = 0;

 /* 
 // no arguments given
 */
 if(argc == 1) {
  //fprintf(stderr, "This program needs arguments....\n\n");
  //print_help(1);
 }

 while((opt = getopt(argc, argv, "hVvgf:pt:m:co:")) != -1) {
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
    case 'f':
      fileName  = strdup(optarg);
      //multMatrices(optarg);
      break;
    case 'c':
      multiply  = 1;
      break;
    case 'p':
      print   = 1;
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

  if (multiply && fileName)
    multMatrices(fileName, nthrds, method, print);  

 return 0;
}
