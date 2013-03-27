/**
 * \file   mat-mult.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   February 2013
 * \brief  Header file for dense matrix multiplication.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "matrix.h"

#ifdef __F4RT_HAVE_PTHREAD_H
#include "pthrd/mat-mult-pthrd.h"
#endif
#ifdef __F4RT_HAVE_INTEL_TBB
#include "tbb/mat-mult-tbb.h"
#endif
#ifdef __F4RT_HAVE_OPENMP
#include "omp/mat-mult-omp.h"
#endif
#if defined(__F4RT_HAVE_KAAPI)
#include "kaapi/mat-mult-kaapi.h"
#endif
#include "seq/mat-mult-seq.h"

void prepareMult(Matrix& A, Matrix& B, char* str);

void multiply(
  Matrix& C, const Matrix& A, const Matrix& B, const int nthrds,
  const int blocksize, const int method, const int dimension, 
  const int affinity, int impose, int outerloop
  );

void multMatrices(
  char* str1, char* str2, int nthrds, int method, int affinity,
  int blocksize, int dimension, int impose, int outerloop, int print
  );

void multEqualMatrices(
  char* str, int nthrds, int method, int affinity, int blocksize, 
  int dimension, int impose, int outerloop, int print
  );
