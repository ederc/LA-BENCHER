/**
 * \file   mat-mult.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   February 2013
 * \brief  General source file for dense matrix multiplication.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-mult.h"

void prepareMult(Matrix& A, Matrix& B, char* str) {
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);

  // let B be just a copy of A
  // we will then multiply A*B^T
  B.transpose(A);
}

void multiply(Matrix& C, const Matrix& A, const Matrix& B, const int nthrds,
              const uint32 blocksize, const int method, const int dimension, 
              const int affinity, int impose, int outerloop) {
  // C = A*B^T
  printf("Impose %d\n",impose);
  if (method == 2) { // TBB
#ifdef __F4RT_HAVE_INTEL_TBB
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
    }
    if (dimension == 2) {
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
    /*
    if (dimension == 3) {
      if (affinity == 1) {
        multTBBAffine3d(C, A, B, nthrds, blocksize, impose);
      } else {
        if (affinity == 2) {
          multTBBSimple3d(C, A, B, nthrds, blocksize, impose);
        } else {
          multTBBAuto3d(C, A, B, nthrds, blocksize, impose);
        }
      }
    }
    */
#else
    multSEQ(C, A, B, blocksize, impose);
#endif
  }
  if (method == 1) { // OpenMP
#ifdef __F4RT_HAVE_OPENMP    
    if (dimension == 1) {
      if (outerloop == 1)
        multOMP1dOuter(C, A, B, nthrds, blocksize, impose);
      else
        multOMP1dInner(C, A, B, nthrds, blocksize, impose);
    }
    if (dimension == 2)
      multOMP2d(C, A, B, nthrds, blocksize, impose);
#else
    multSEQ(C, A, B, blocksize, impose);
#endif
  }
  if (method == 0) // plain sequential w/o scheduler overhead
    multSEQ(C, A, B, blocksize, impose);
  if (method == 3) {// xkaapi 
#if defined(__F4RT_HAVE_KAAPIC)
    // TODO: How to enlarge blocksize without corrupting the computation?
    if (dimension == 1) {
      multKAAPIC1d(C, A, B, nthrds, blocksize, impose);
    }
    if (dimension == 2) {
      multKAAPIC2d(C, A, B, nthrds, blocksize, impose);
    }
#else
    multSEQ(C, A, B, blocksize, impose);
#endif
  }
  if (method == 4) { // pthreads
#ifdef __F4RT_HAVE_PTHREAD_H
    multPT(C, A, B, nthrds, blocksize, impose);
#else
    multSEQ(C, A, B, blocksize, impose);
#endif
  }
}

void multMatrices(char* str1, char* str2, int nthrds, int method, int affinity, 
                  uint32 blocksize, int dimension, int impose, int outerloop, 
                  int print) {
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

  multiply(C, A, B, nthrds, blocksize, method, dimension, affinity, impose, outerloop);

  if (print)
    C.print();
  // clear memory
  A.clear();
  B.clear();
  C.clear();
}

void multEqualMatrices( char* str, int nthrds, int method, int affinity,
                        uint32 blocksize, int dimension, int impose, 
                        int outerloop, int print) {
  Matrix A, B;

  // read files, stores matrices, etc
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);

  // let B be just a copy of A
  // we will then multiply A*B^T
  B.transpose(A);
  
  Matrix C(A.nRows(), B.nRows());

  multiply(C, A, B, nthrds, blocksize, method, dimension, affinity, impose, outerloop);

  if (print)
    C.print();
  // clear memory
  A.clear();
  B.clear();
  C.clear();
}
