/**
 * \file   mat-mult-seq.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for sequential dense matrix multiplication.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_MULT_SEQ_H
#define F4RT_MAT_MULT_SEQ_H

#include <matrix.h>

// multiplies A*B^T and stores it in *this
void multSEQ(Matrix& C, const Matrix& A, const Matrix& B, int blocksize, int impose);
#endif
