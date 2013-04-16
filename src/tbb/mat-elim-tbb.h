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

#ifdef __F4RT_HAVE_INTEL_TBB
void elimTBB(Matrix& A, uint32 blocksize);

void elimNaiveTBBModP1dAuto(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);
void elimNaiveTBBModP1dAffine(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);
void elimNaiveTBBModP1dSimple(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

void elimNaiveTBBModP2dAuto(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);
void elimNaiveTBBModP2dAffine(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);
void elimNaiveTBBModP2dSimple(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

void elimNaiveTBBModP1dAutoPivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);
void elimNaiveTBBModP1dAffinePivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);
void elimNaiveTBBModP1dSimplePivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

void elimNaiveTBBModP2dAutoPivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);
void elimNaiveTBBModP2dAffinePivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);
void elimNaiveTBBModP2dSimplePivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

// cache-oblivious stuff

void elimCoTBBBaseModP( mat *M, const uint32 k1, const uint32 i1,
                        const uint32 j1, const uint32 rows, const uint32 cols,
                        uint64 size, uint64 prime, mat *neg_inv_piv,
                        int nthrds, uint32 blocksize);

void ATBB( mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void B1TBB(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void C1TBB(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void D1TBB(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void elimCoTBBModP(Matrix& M, int nthrds, uint32 blocksize, uint64 prime);
#endif
#endif
