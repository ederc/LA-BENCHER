/**
 * \file   mat-elim.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_H
#define F4RT_MAT_ELIM_H

#include "matrix.h"
#include "seq/mat-elim-seq.h"

#ifdef __F4RT_HAVE_PTHREAD_H
#include "pthrd/mat-elim-pthrd.h"
#endif
#ifdef __F4RT_HAVE_INTEL_TBB
#include "tbb/mat-elim-tbb.h"
#endif
#ifdef __F4RT_HAVE_OPENMP
#include "omp/mat-elim-omp.h"
#endif
#ifdef __F4RT_HAVE_KAAPIC
#include "kaapi/mat-elim-kaapi.h"
#endif
#ifdef __F4RT_HAVE_LAPACK
#include "blas/mat-elim-blas.h"
#endif


void eliminate(Matrix& A, const int nthrds, const uint32 blocksize, 
              const int method, const int dimension, const int affinity, 
              int outerloop, int pivoting, int cacheOblivious, uint64 prime);

void eliminateMatrix( char* str, int nthrds, int method, int affinity, 
                      uint32 blocksize, int dimension, int outerloop, 
                      int pivoting, int cacheOblivious, uint64 prime, int print);
#endif
