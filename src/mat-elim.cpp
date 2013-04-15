/**
 * \file   mat-elim.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  General source file for Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim.h"

void eliminate(Matrix& A, const int nthrds, const uint32 blocksize, 
              const int method, const int dimension, const int affinity, 
              int outerloop, int pivoting, int cacheOblivious, uint64 prime) {
  // sequential method
  if (method == 0) {
    if (cacheOblivious == 0) {
      if (pivoting == 1) {
        elimNaiveSEQModPPivot(A, blocksize, prime);
      } else {
        elimNaiveSEQModP(A, blocksize, prime);
      }
      //A.print();
    } else {
      elimCoSEQModP(A, blocksize, prime);
    }
  }
  // OpenMP
  if (method == 1) {
#ifdef __F4RT_HAVE_OPENMP    
    if (dimension == 1) {
      if (cacheOblivious == 0) {
        if (pivoting == 1) {
          elimNaiveOMPModP1dOuterPivot(A, nthrds, blocksize, prime);
        } else {
          elimNaiveOMPModP1dOuter(A, nthrds, blocksize, prime);
        }
      } else {
        elimCoOMPModP(A, nthrds, blocksize, prime);
      }
    }
#else
    elimNaiveSEQModP(A, blocksize, prime);
#endif
  }
  // TBB
  if (method == 2) {
#ifdef __F4RT_HAVE_INTEL_TBB
    if (cacheOblivious == 0) {
      if (dimension == 1) {
        if (affinity == 1) {
          if (pivoting == 1) {
            elimNaiveTBBModP1dAffinePivot(A, nthrds, blocksize, prime);
          } else {
            elimNaiveTBBModP1dAffine(A, nthrds, blocksize, prime);
          }
        } else {
          if (affinity == 2) {
            if (pivoting == 1) {
              elimNaiveTBBModP1dSimplePivot(A, nthrds, blocksize, prime);
            } else {
              elimNaiveTBBModP1dSimple(A, nthrds, blocksize, prime);
            }
          } else {
            if (pivoting == 1) {
              elimNaiveTBBModP1dAutoPivot(A, nthrds, blocksize, prime);
            } else {
              elimNaiveTBBModP1dAuto(A, nthrds, blocksize, prime);
            }
          }
        }
      }
      if (dimension == 2) {
        if (affinity == 1) {
          if (pivoting == 1) {
            elimNaiveTBBModP2dAffinePivot(A, nthrds, blocksize, prime);
          } else {
            elimNaiveTBBModP2dAffine(A, nthrds, blocksize, prime);
          }
        } else {
          if (affinity == 2) {
            if (pivoting == 1) {
              elimNaiveTBBModP2dSimplePivot(A, nthrds, blocksize, prime);
            } else {
              elimNaiveTBBModP2dSimple(A, nthrds, blocksize, prime);
            }
          } else {
            if (pivoting == 1) {
              elimNaiveTBBModP2dAutoPivot(A, nthrds, blocksize, prime);
            } else {
              elimNaiveTBBModP2dAuto(A, nthrds, blocksize, prime);
            }
          }
        }
      }
    } else {
      if (dimension == 1) {
        if (affinity == 1) {
          // TODO
        } else {
          if (affinity == 2) {
            // TODO
          } else {
            elimCoTBBModP(A, nthrds, blocksize, prime);
          }
        }
      }
    }
#else
    elimNaiveSEQModP(A, blocksize, prime);
#endif
  }
  // KAAPI
  if (method == 3) {
#ifdef __F4RT_HAVE_KAAPI
    if (pivoting == 1) {
      elimNaiveKAAPICModP1dPivot(A, nthrds, blocksize, prime);
    } else {
      elimNaiveKAAPICModP1d(A, nthrds, blocksize, prime);
    }
#else
    elimNaiveSEQModP(A, blocksize, prime);
#endif
  }
  // pthread
  if (method == 4) {
#ifdef __F4RT_HAVE_PTHREAD_H
    if (pivoting == 1) {
      elimNaivePTHRDModP1dPivot(A, nthrds, blocksize, prime);
    } else {
      elimNaivePTHRDModP1d(A, nthrds, blocksize, prime);
    }
#else
    elimNaiveSEQModP(A, blocksize, prime);
#endif
  }
  // only for testing purposes, no user option
  if (method == 5) {
      elimCoSEQModPOld(A, blocksize, prime);
  }
}

void eliminateMatrix( char* str, int nthrds, int method, int affinity,
                      uint32 blocksize, int dimension, int outerloop, 
                      int pivoting, int cacheOblivious, uint64 prime, 
                      int print) {
  Matrix A;

  // read files, stores matrices, etc
  FILE* file  = fopen(str,"rb");
  // take A from file
  A.read(file);
 
  eliminate(A, nthrds, blocksize, method, dimension, affinity,
            outerloop, pivoting, cacheOblivious, prime);

  if (print)
    A.print();
  // clear memory
  A.clear();
}
