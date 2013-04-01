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

mat negInverseModP(mat a, uint64 prime);

void cleanUpModP(Matrix& A, uint64 prime);

void elimSEQ(Matrix& A, int blocksize);
void elimNaiveSEQModP(Matrix& A, int blocksize, uint64 prime);
#endif
