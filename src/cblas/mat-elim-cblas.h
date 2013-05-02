/**
 * \file   mat-elim-seq.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for sequential Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_SEQ_H
#define F4RT_MAT_ELIM_SEQ_H

#include <matrix.h>
#include "../mat-elim-tools.h"


void elimNaiveSEQModP(Matrix& A, uint32 blocksize, uint64 prime);
void elimNaiveSEQModPPivot(Matrix& A, uint32 blocksize, uint64 prime);

void A( mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime,
        mat *neg_inv_piv, uint32 blocksize);

void B1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime,
        mat *neg_inv_piv, uint32 blocksize);

void C1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime,
        mat *neg_inv_piv, uint32 blocksize);

void D1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime,
        mat *neg_inv_piv, uint32 blocksize);

void elimCoSEQBaseModP( mat *M, const uint32 k1, const uint32 i1, 
                        const uint32 j1, const uint32 rows, const uint32 cols,
                        uint64 size, uint64 prime,
                        mat *neg_inv_piv, uint32 blocksize);
void elimCoSEQModPOld(Matrix& A, uint32 blocksize, uint64 prime);
void elimCoSEQModP(Matrix& A, uint32 blocksize, uint64 prime);
#endif
