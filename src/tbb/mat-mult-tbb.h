/**
 * \file   mat-mult-tbb.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for dense matrix multiplication using Intel TBB.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_MULT_TBB_H
#define F4RT_MAT_MULTTBB_H

#include <tbb/tbb.h>
#include <matrix.h>

// multiplies A*B^T and stores it in *this
void multTBBAuto( Matrix& C, const Matrix& A, const Matrix& B, int nthrds,
                  int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multTBBAffine( Matrix& C, const Matrix& A, const Matrix& B, int nthrds,
                    int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multTBBSimple( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                    int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multTBBAuto2d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                    int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multTBBAffine2d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                      int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multTBBSimple2d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                      int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multTBBAuto3d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                    int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multTBBAffine3d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                      int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multTBBSimple3d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                      int blocksize, int impose);
#endif
