#ifndef F4RT_MATRIX_OMP_H
#define F4RT_MATRIX_OMP_H

#include <matrix.h>

// multiplies A*B^T and stores it in *this
void multOMP(Matrix& C, const Matrix& A, const Matrix& B, int nthrds);
#endif
