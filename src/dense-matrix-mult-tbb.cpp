#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "matrix-tbb.h"

#define PACKAGE "dense-mult-tbb"
#define VERSION "0.0.1"


void print_help(int exval) {
 printf("%s,%s show working getopt example\n", PACKAGE, VERSION); 
 printf("%s [-h] [-V] [-f FILE] [-o FILE]\n\n", PACKAGE);

 printf("  -h              print this help and exit\n");
 printf("  -V              print version and exit\n\n");

 printf("  -v              set verbose flag\n");
 printf("  -g              generate a new random uint16 matrix\n");
 printf("  -f FILE         set intput file\n");
 printf("  -m              if input file is set, multiply matrix with its own transpose\n");
 printf("  -p              if matrix multiplication took place, print of resulting matrix\n");
 printf("                  (no printing of resulting matrix by default)\n");

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
  std::cout << str << std::endl;
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);

  // let B be just a copy of A
  // we will then multiply A*B^T
  B.copy(A);
}

void multMatrices(char* str, int print) {
  Matrix A, B;

  // read files, stores matrices, etc
  prepareMult(A, B, str);
  
  Matrix C(A.nRows(), B.nRows());

  // C = A*B^T
  multTBB(C, A, B);
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
 int print = 0, multiply  = 0;;

 /* 
 // no arguments given
 */
 if(argc == 1) {
  //fprintf(stderr, "This program needs arguments....\n\n");
  //print_help(1);
 }

 while((opt = getopt(argc, argv, "hVvgf:pmo:")) != -1) {
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
    case 'm':
      multiply  = 1;
      break;
    case 'p':
      print   = 1;
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
    multMatrices(fileName, print);  

 return 0;
}
/*
# include "stdio.h"
# include "stdint.h"
main()
{
    uint16_t m1[10000][10],i,j,k,m2[10][10],mult[10][10],r1,c1,r2,c2;
    printf("Enter number of rows and columns of first matrix (less than 10)\n");
    scanf("%d%d",&r1,&c1);
    printf("Enter number of rows and columns of second matrix (less than 10)\n");
    scanf("%d%d",&r2,&c2);
    if(r2==c1)
    {
        printf("Enter rows and columns of First matrix \n");
        printf("Row wise\n");
        for(i=0;i<r1;i++)
            for(j=0;j<c1;j++)
                scanf("%d",&m1[i][j]);
        printf("First Matrix is :\n");
        for(i=0;i<r1;i++)
        {
            for(j=0;j<c1;j++)
                printf("%d\t",m1[i][j]);
            printf("\n");
        }
        printf("Enter rows and columns of Second matrix \n");
        printf("Row wise\n");
        for(i=0;i<r2;i++)
            for(j=0;j<c2;j++)
                scanf("%d",&m2[i][j]);
        printf("Second Matrix is:\n");
        for(i=0;i<r2;i++)
        {
            for(j=0;j<c2;j++)
                printf("%d\t",m2[i][j]);
            printf("\n");
        }
        printf("Multiplication of the Matrices:\n");
        for(i=0;i<r1;i++)
        {
            for(j=0;j<c2;j++)
            {
                mult[i][j]=0;
                for(k=0;k<r1;k++)
                    mult[i][j]+=m1[i][k]*m2[k][j];
                printf("%d\t",mult[i][j]);
            }
            printf("\n");
        }
    }
    else
    {
        printf("Matrix multiplication cannot be done");
    }
    return 0;
}
*/
