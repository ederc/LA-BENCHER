#ifndef F4RT_MATRIX_KAAPI_H
#define F4RT_MATRIX_KAAPI_H

#include <matrix.h>

#if defined(__F4RT_HAVE_KAAPI) && defined(__F4RT_ENABLE_KAAPI)
static int BLOCSIZE = 1;

/********************************
 * 1 Dimensional implementation
 *******************************/

// Parallel bloc matrix product
struct TaskMatProduct: public ka::Task<0>::Signature{};

template<>
struct TaskBodyCPU<TaskMatProduct> {
  void operator()( 
    Matrix& C,
    const Matrix& A,
    const Matrix& B,
    const int nthrds,
    const int blocksize,
    const int impose
    ) {
    uint32 l, m, n;
    if (impose == 1) {
      l = A.nRows();
      m = B.nRows();
      n = B.nCols();
    } else {
      l = A.nRows();
      m = B.nCols();
      n = B.nRows();
    }
    std::cout << "Matrix Multiplication" << std::endl;
#if __F4RT_DEBUG
    std::cout << std::endl;
    std::cout << "A => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
    std::cout << "B => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
    std::cout << "C => " << C.nRows() << "-" << C.nCols() << "-" << C.nEntries() << std::endl;
#endif
    if (impose == 1) {
      for (size_t i=0; i<l; i += blocksize) {
        ka::rangeindex ri(i, i+blocksize);
        for (size_t j=0; j<m; j += blocksize) {
          ka::rangeindex rj(j, j+blocksize);
          float sum = 0;
          for (size_t k=0; k<n; k += 1) {
            sum += A.entries[k+i*n] * B.entries[k+j*n];
          }
          C.entries[j+i*m]  = (float) (sum);
        }
      }
    } else {
      for (size_t i=0; i<l; i += blocksize) {
        ka::rangeindex ri(i, i+blocksize);
        for (size_t j=0; j<m; j += blocksize) {
          ka::rangeindex rj(j, j+blocksize);
          float sum = 0;
          for (size_t k=0; k<n; k += 1) {
            sum += A.entries[k+i*n] * B.entries[j+k*m];
          }
          C.entries[j+i*m]  = (float) (sum);
        }
      }
    }
  }
};

/********************************
 * 2 Dimensional implementation
 *******************************/
/*
// Parallel bloc matrix product (lower level)
struct TaskMatProductL2: public ka::Task<3>::Signature<
      ka::R<ka::range2d<double> >,
      ka::R<ka::range2d<double> >,
      ka::RPWP<ka::range2d<double> >
>{};

template<>
struct TaskBodyCPU<TaskMatProductL2> {
  void operator()(
    ka::range2d_r<double> A,
    ka::range2d_r<double> B,
    ka::range2d_rpwp<double> C
  ) {
    size_t M = A->dim(0);
    size_t K = B->dim(0);
    size_t N = B->dim(1);
    int bloc = BLOCSIZE/4;

    for (size_t i=0; i<M; i += bloc) {
      ka::rangeindex ri(i, i+bloc);
      for (size_t j=0; j<N; j += bloc) {
        ka::rangeindex rj(j, j+bloc);
        for (size_t k=0; k<K; k += bloc) {
          ka::rangeindex rk(k, k+bloc);
          ka::Spawn<TaskDGEMM>() (
            CblasNoTrans, CblasNoTrans, 1.0, A(ri,rk), B(rk,rj), 1.0, C(ri,rj)
          );
        }
      }
    }
  }
};

// Parallel bloc matrix product (upper level)
struct TaskMatProductL1: public ka::Task<3>::Signature<
      ka::R<ka::range2d<double> >,
      ka::R<ka::range2d<double> >,
      ka::RPWP<ka::range2d<double> >
>{};

template<>
struct TaskBodyCPU<TaskMatProductL1> {
  void operator()(
    ka::range2d_r<double> A,
    ka::range2d_r<double> B,
    ka::range2d_rpwp<double> C
  ) {
    size_t M = A->dim(0);
    size_t K = B->dim(0);
    size_t N = B->dim(1);
    int bloc = BLOCSIZE;

    for (size_t i=0; i<M; i += bloc) {
      ka::rangeindex ri(i, i+bloc);
      for (size_t j=0; j<N; j += bloc) {
        ka::rangeindex rj(j, j+bloc);
        for (size_t k=0; k<K; k += bloc) {
          ka::rangeindex rk(k, k+bloc);
          ka::Spawn<TaskMatProduct2>(ka::SetStaticSched())(
            A(ri,rk), B(rk,rj), C(ri,rj) 
          );
        }
      }
    }
  }
};
*/
// multiplies A*B^T and stores it in *this
void multKAAPI(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose);

#endif
#endif
