#ifndef F4RT_OMP_H
#define F4RT_OMP_H

#include <matrix.h>

// multiplies A*B^T and stores it in *this
void multOMP1dInner(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multOMP1dOuter(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multOMP2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose);
#endif
