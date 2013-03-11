#ifndef F4RT_MATRIX_SEQ_H
#define F4RT_MATRIX_SEQ_H

#include <matrix.h>

// multiplies A*B^T and stores it in *this
void multSEQ(Matrix& C, const Matrix& A, const Matrix& B, int blocksize, int impose);
#endif
