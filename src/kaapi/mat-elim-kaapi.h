/**
 * \file   mat-elim-kaapi.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   April 2013
 * \brief  Header file for Gaussian Elimination using XKAAPI.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef F4RT_MAT_ELIM_KAAPI_H
#define F4RT_MAT_ELIM_KAAPI_H

//#include <matrix.h>
#include "../mat-elim-tools.h"

#if defined(__F4RT_HAVE_KAAPIC)

void elimKAAPIC(mat *a_entries, uint32 rows, uint32 cols, uint32 blocksize);

void elimNaiveKAAPICModP1d(mat *a_entries, uint32 rows, uint32 cols, int nthrds, uint32 blocksize, uint64 prime);

void elimNaiveKAAPICModP1dPivot(mat *a_entries, uint32 rows, uint32 cols, int nthrds, uint32 blocksize, uint64 prime);

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

void elimCoKAAPICModP(mat *a_entries, uint32 rows, uint32 cols, int nthrds, uint32 blocksize, uint64 prime);
#endif
#endif
#ifdef __cplusplus
}
#endif
