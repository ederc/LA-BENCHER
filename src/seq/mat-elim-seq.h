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

void elimNaiveSEQModP(Matrix& A, int blocksize, uint64 prime);
void elimNaiveSEQModPPivot(Matrix& A, int blocksize, uint64 prime);

void A( mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime, mat *neg_inv_piv);

void B1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime, mat *neg_inv_piv);

void B2(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime, mat *neg_inv_piv);

void C1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime, mat *neg_inv_piv);

void C2(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime, mat *neg_inv_piv);

void D1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime, mat *neg_inv_piv);

void elimCoSEQModP(Matrix& A, int blocksize, uint64 prime);
void elimCoSEQModPReordered(Matrix& A, int blocksize, uint64 prime);
#endif
