/**
 * \file   mat-elim-tbb.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for Gaussian Elimination using Intel TBB.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_TBB_H
#define F4RT_MAT_ELIM_TBB_H

#include <matrix.h>
#include "../mat-elim-tools.h"

void elimTBB(Matrix& A, int blocksize);

void elimNaiveTBBModP1dAuto(Matrix& A, int nthrds, int blocksize, uint64 prime);
void elimNaiveTBBModP1dAffine(Matrix& A, int nthrds, int blocksize, uint64 prime);
void elimNaiveTBBModP1dSimple(Matrix& A, int nthrds, int blocksize, uint64 prime);

void elimNaiveTBBModP2dAuto(Matrix& A, int nthrds, int blocksize, uint64 prime);
void elimNaiveTBBModP2dAffine(Matrix& A, int nthrds, int blocksize, uint64 prime);
void elimNaiveTBBModP2dSimple(Matrix& A, int nthrds, int blocksize, uint64 prime);

void elimNaiveTBBModP1dAutoPivot(Matrix& A, int nthrds, int blocksize, uint64 prime);
void elimNaiveTBBModP1dAffinePivot(Matrix& A, int nthrds, int blocksize, uint64 prime);
void elimNaiveTBBModP1dSimplePivot(Matrix& A, int nthrds, int blocksize, uint64 prime);

void elimNaiveTBBModP2dAutoPivot(Matrix& A, int nthrds, int blocksize, uint64 prime);
void elimNaiveTBBModP2dAffinePivot(Matrix& A, int nthrds, int blocksize, uint64 prime);
void elimNaiveTBBModP2dSimplePivot(Matrix& A, int nthrds, int blocksize, uint64 prime);
#endif
