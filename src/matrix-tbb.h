#ifndef F4RT_MATRIX_TBB_H
#define F4RT_MATRIX_TBB_H

#include <tbb/tbb.h>
#include <matrix.h>

// multiplies A*B^T and stores it in *this
void multTBBAuto(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize);

// multiplies A*B^T and stores it in *this
void multTBBAffine(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize);

// multiplies A*B^T and stores it in *this
void multTBBAuto2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize);

// multiplies A*B^T and stores it in *this
void multTBBAffine2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize);
#endif
