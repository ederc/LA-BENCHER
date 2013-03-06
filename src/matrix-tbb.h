#ifndef F4RT_MATRIX_TBB_H
#define F4RT_MATRIX_TBB_H

#include <tbb/tbb.h>
#include <matrix.h>
//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range2d.h>

using namespace tbb;

class MatrixMultiply2D {
  Matrix _a, _b, _c;
  uint32 _a_rows, _b_rows, _cols;
  public:
  void operator()( const blocked_range2d<size_t>& r ){
    //Matrix __a = _a, __b = _b, __c = _c;
    for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i )
      for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ) {
        uint32 sum = 0;
        for( size_t k=0; k<_cols; ++k )
          sum += _a(i,k) * _b(j,k);
        _c(i,j) = sum;
      }
  }
  MatrixMultiply2D(Matrix& C, const Matrix& A, const Matrix& B) :
    _a(A), _b(B), _c(C), _a_rows(A.nRows()), _b_rows(B.nRows()), _cols(B.nCols()) {}
};

class MatrixMultiply1D {
  Matrix _a, _b, _c;
  uint32 _a_rows, _b_rows, _cols;
  public:
  void operator()( const blocked_range<size_t>& r ) {
    //Matrix A = _a, B = _b, C = _c;
    for( size_t i=r.begin(); i!=r.end(); ++i )
      for( size_t j=_b_rows; j!=_b_rows; ++j ) {
        uint32 sum = 0;
        for( size_t k=0; k<_cols; ++k )
          sum += _a(i,k) * _b(j,k);
        _c(i,j) = sum;
      }
  }
  MatrixMultiply1D(Matrix& C, const Matrix& A, const Matrix& B) :
    _a(A), _b(B), _c(C), _a_rows(A.nRows()), _b_rows(B.nRows()), _cols(B.nCols()) {}
};

// multiplies A*B^T and stores it in *this
void multTBB(Matrix& C, const Matrix& A, const Matrix& B, int nthrds);
#endif
