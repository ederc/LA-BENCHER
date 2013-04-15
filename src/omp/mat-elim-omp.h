/**
 * \file   mat-elim-omp.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for Gaussian Elimination using OMP.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_OMP_H
#define F4RT_MAT_ELIM_OMP_H

#include <matrix.h>
#include "../mat-elim-tools.h"

void elimOMP(Matrix& A, uint32 blocksize);
void elimNaiveOMPModP1dOuter(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

void elimNaiveOMPModP1dOuterPivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

// cache-oblivious stuff

void elimCoOMPBaseModP( mat *M, const uint32 k1, const uint32 i1,
                        const uint32 j1, const uint32 rows, const uint32 cols,
                        uint64 size, uint64 prime, mat *neg_inv_piv,
                        int nthrds, uint32 blocksize);

void AOMP( mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void B1OMP(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void C1OMP(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void D1OMP(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void elimCoOMPModP(Matrix& M, int nthrds, uint32 blocksize, uint64 prime);
#endif
