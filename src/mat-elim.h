/**
 * \file   mat-elim.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   February 2013
 * \brief  Header file for Gaussian Elimination.
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
#include "seq/mat-elim-seq.h"

void eliminateMatrix( 
  char* str, int nthrds, int method, int affinity, 
  int blocksize, int dimension, int outerloop, int print
  );
