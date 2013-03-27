/**
 * \file   mat-elim.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   February 2013
 * \brief  General source file for Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim.h"

#ifdef __F4RT_HAVE_PTHREAD_H
#include "pthrd/mat-elim-pthrd.h"
#endif
#ifdef __F4RT_HAVE_INTEL_TBB
#include "tbb/mat-elim-tbb.h"
#endif
#ifdef __F4RT_HAVE_OPENMP
#include "omp/mat-elim-omp.h"
#endif
#if defined(__F4RT_HAVE_KAAPI)
#include "kaapi/mat-elim-kaapi.h"
#endif
#include "seq/mat-elim-seq.h"
