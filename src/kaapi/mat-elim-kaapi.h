/**
 * \file   mat-elim-kaapi.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   April 2013
 * \brief  Header file for Gaussian Elimination using XKAAPI.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_KAAPI_H
#define F4RT_MAT_ELIM_KAAPI_H

#include <matrix.h>
#include "../mat-elim-tools.h"

#if defined(__F4RT_HAVE_KAAPI)

void elimKAAPIC(Matrix& A, uint32 blocksize);

void elimNaiveKAAPICModP1d(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

void elimNaiveKAAPICModP1dPivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime);

// cache-oblivious stuff

void elimCoKAAPICBaseModP( mat *M, const uint32 k1, const uint32 i1,
                        const uint32 j1, const uint32 rows, const uint32 cols,
                        uint64 size, uint64 prime, mat *neg_inv_piv,
                        int nthrds, uint32 blocksize);

void AKAAPIC( mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void B1KAAPIC(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void C1KAAPIC(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void D1KAAPIC(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize);

void elimCoKAAPICModP(Matrix& M, int nthrds, uint32 blocksize, uint64 prime);
#endif
#endif
