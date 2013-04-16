/**
 * \file   mat-elim-pthrd.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   April 2013
 * \brief  Header file for Gaussian Elimination using pThread.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_PTHRD_H
#define F4RT_MAT_ELIM_PTHRD_H

#include <matrix.h>
#include "../mat-elim-tools.h"

#ifdef __F4RT_HAVE_PTHREAD_H

struct paramsElim {
  mat *a;
  int size;
  mat inv;
  uint64 prime;
  uint32 index;
  uint32 start;
  uint32 m;
  uint32 n;
  int tid;
};


void elimPTHRD(Matrix& A, uint32 blocksize);

void elimNaivePTHRDModP1d(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

void elimNaivePTHRDModP1dPivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

#endif
#endif
