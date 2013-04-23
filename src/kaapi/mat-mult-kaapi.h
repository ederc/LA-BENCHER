/**
 * \file   mat-mult-kaapi.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for dense matrix multiplication using KAAPI.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_MULT_KAAPI_H
#define F4RT_MAT_MULT_KAAPI_H

#include <matrix.h>

#if defined(__F4RT_HAVE_KAAPIC)

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
void multKAAPIC1d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, uint32 blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multKAAPIC2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, uint32 blocksize, int impose);

// multiplies A*B^T and stores it in *this
void multKAAPI(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, uint32 blocksize, int impose);

#endif
#endif
