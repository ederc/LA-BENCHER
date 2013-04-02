/**
 * \file   mat-elim.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  General source file for Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim.h"

void eliminate(Matrix& A, const int nthrds, const int blocksize, 
              const int method, const int dimension, const int affinity, 
              int outerloop, uint64 prime) {
  // sequential method
  if (method == 0) {
    elimNaiveSEQModP(A, blocksize, prime);
    //A.print();
  }
  // OpenMP
  if (method == 1) {
#ifdef __F4RT_HAVE_OPENMP    
    if (dimension == 1) {
      elimNaiveOMPModP1dOuter(A, nthrds, blocksize, prime);
    }
#else
    elimNaiveSEQModP(A, blocksize, prime);
#endif
  }
  // TBB
  if (method == 2) {
#ifdef __F4RT_HAVE_INTEL_TBB
    if (dimension == 1) {
      if (affinity == 1) {
        elimNaiveTBBModP1dAffine(A, nthrds, blocksize, prime);
      } else {
        if (affinity == 2) {
          elimNaiveTBBModP1dSimple(A, nthrds, blocksize, prime);
        } else {
          elimNaiveTBBModP1dAuto(A, nthrds, blocksize, prime);
        }
      }
    }
    if (dimension == 2) {
      if (affinity == 1) {
        elimNaiveTBBModP2dAffine(A, nthrds, blocksize, prime);
      } else {
        if (affinity == 2) {
          elimNaiveTBBModP2dSimple(A, nthrds, blocksize, prime);
        } else {
          elimNaiveTBBModP2dAuto(A, nthrds, blocksize, prime);
        }
      }
    }
#else
    elimNaiveSEQModP(A, blocksize, prime);
#endif
  }
  // KAAPI
  if (method == 3) {
  }
  // pthread
  if (method == 4) {
  }
}

void eliminateMatrix( char* str, int nthrds, int method, int affinity,
                      int blocksize, int dimension, int outerloop, 
                      uint64 prime, int print) {
  Matrix A;

  // read files, stores matrices, etc
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);
 
  eliminate(A, nthrds, blocksize, method, dimension, affinity, outerloop, prime);

  if (print)
    A.print();
  // clear memory
  A.clear();
}
