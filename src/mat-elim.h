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

#ifdef __F4RT_HAVE_PTHREAD_H
#include "pthrd/mat-mult-pthrd.h"
#endif
#ifdef __F4RT_HAVE_INTEL_TBB
#include "tbb/mat-elim-tbb.h"
#endif
#ifdef __F4RT_HAVE_OPENMP
#include "omp/mat-elim-omp.h"
#endif
#if defined(__F4RT_HAVE_KAAPI)
#include "kaapi/mat-mult-kaapi.h"
#endif
#include "seq/mat-elim-seq.h"

void eliminate(Matrix& A, const int nthrds, const int blocksize, 
              const int method, const int dimension, const int affinity, 
              int outerloop, uint64 prime);

void eliminateMatrix( char* str, int nthrds, int method, int affinity, 
                      int blocksize, int dimension, int outerloop, 
                      uint64 prime, int print);
#endif
