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
  }
  // TBB
  if (method == 2) {
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
