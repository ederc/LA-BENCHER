/**
 * \file   mat-mul-omp.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for dense matrix multiplication using OpenMP.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_MULT_OMP_H
#define F4RT_MAT_MULT_OMP_H

#include <matrix.h>

// multiplies A*B^T and stores it in *this
void multOMP1dInner(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multOMP1dOuter(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multOMP2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose);
#endif
