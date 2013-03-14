#ifndef F4RT_MATRIX_KAAPI_H
#define F4RT_MATRIX_KAAPI_H

#include <matrix.h>

#if defined(__F4RT_HAVE_KAAPI) && defined(__F4RT_ENABLE_KAAPI)
static int BLOCSIZE = 0;


/********************************************************
 * start: from matrix.h in examples of xkaapi repository
 *******************************************************/
/* Note: this file defines 
    - tasks for some BLAS-3 over double/float 
    - tasks for some LAPACK functions.
  The LAPACK C interface is assumed to be the ATLAS version
  where it exist both version of LU or LLt factorization for
  row major or column major representation of matrix.
  
  Here, we assume row major representation of matrix. The
*/


/* Task Print Matrix
 * this task prints the matrix using Maple matrix format
 */
template<typename T>
struct TaskPrintMatrix : public ka::Task<2>::Signature<
  std::string, 
  ka::R<ka::range2d<T> 
> > {};

template<typename T>
struct TaskBodyCPU<TaskPrintMatrix<T> > {
  void operator() ( std::string msg, ka::range2d_r<T> A  )
  {
    size_t d0 = A->dim(0);
    size_t d1 = A->dim(1);
    std::cout << msg << " :=matrix( [" << std::endl;
    for (size_t i=0; i < d0; ++i)
    {
      std::cout << "[";
      for (size_t j=0; j < d1; ++j)
      {
        std::cout << std::setw(18) << std::setprecision(15) << std::scientific << A(i,j) << (j == d1-1 ? "" : ", ");
      }
      std::cout << "]" << (i == d0-1 ? ' ' : ',') << std::endl;
    }
    std::cout << "]);" << std::endl;
  }
};


/* =================== CBLAS routines =================== */

/* Solve : L(rk,rk) * X =  * A(rk,rj) 
    ie:  X <- L(rk,rk)^-1 * A(rk,rj) 
*/
struct TaskDTRSM_left: public ka::Task<2>::Signature<
      ka::R<ka::range2d<double> >, /* Akk */
      ka::RW<ka::range2d<double> > /* Akj */
>{};
template<>
struct TaskBodyCPU<TaskDTRSM_left> {
  void operator()( ka::range2d_r<double> Akk, ka::range2d_rw<double> Akj )
  {
    const double* const a = Akk->ptr();
    const int lda = Akk->lda();

    double* const b = Akj->ptr();
    const int ldb   = Akj->lda();
    const int n     = Akj->dim(0);
    const int m     = Akj->dim(1);

    cblas_dtrsm
    (
     CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
     n, m, 1., a, lda, b, ldb
    );
  }
};


/* Solve : X * U(rk,rk) =  A(ri,rk) 
    ie:  X <- A(ri,rk) * U(rk,rk)^-1
*/
struct TaskDTRSM_right: public ka::Task<2>::Signature<
      ka::R<ka::range2d<double> >, /* Akk */
      ka::RW<ka::range2d<double> > /* Aik */
>{};
template<>
struct TaskBodyCPU<TaskDTRSM_right> {
  void operator()( ka::range2d_r<double> Akk, ka::range2d_rw<double> Aik )
  {
    const double* const a = Akk->ptr();
    const int lda = Akk->lda();

    double* const b = Aik->ptr();
    const int ldb = Aik->lda();
    const int n = Aik->dim(0); // b.rows();
    const int m = Aik->dim(1); // b.cols();

    cblas_dtrsm
    (
      CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
      n, m, 1., a, lda, b, ldb
    );
  }
};


/* DGEMM rountine to compute
    Aij <- alpha* Aik * Akj + beta Aij
*/
struct TaskDGEMM: public ka::Task<7>::Signature<
      CBLAS_TRANSPOSE,             /* NoTrans/Trans for A */
      CBLAS_TRANSPOSE,             /* NoTrans/Trans for B */
      double,                      /* alpha */
      ka::R<ka::range2d<double> >, /* Aik   */
      ka::R<ka::range2d<double> >, /* Akj   */
      double,                      /* beta */
      ka::RW<ka::range2d<double> > /* Aij   */
>{};
template<>
struct TaskBodyCPU<TaskDGEMM> {
  void operator()
  (
    CBLAS_TRANSPOSE transA,
    CBLAS_TRANSPOSE transB,
    double alpha,
    ka::range2d_r<double> Aik,
    ka::range2d_r<double> Akj,
    double beta,
    ka::range2d_rw<double> Aij
  )
  {
    const double* const a = Aik->ptr();
    const double* const b = Akj->ptr();
    double* const c       = Aij->ptr();

    const int m = Aij->dim(0); 
    const int n = Aij->dim(1);
    const int k = (transA == CblasNoTrans ? Aik->dim(1) : Aik->dim(0) );

    const int lda = Aik->lda();
    const int ldb = Akj->lda();
    const int ldc = Aij->lda();

    cblas_dgemm
    (
      CblasRowMajor, transA, transB,
      m, n, k, alpha, a, lda, b, ldb, beta, c, ldc
    );
  }
};



/* Rank k update
*/
struct TaskDSYRK: public ka::Task<6>::Signature<
  CBLAS_UPLO,                  /* CBLAS Upper / Lower */
  CBLAS_TRANSPOSE,             /* transpose flag */
  double,                      /* alpha */
  ka::R<ka::range2d<double> >, /* A */
  double,                      /* beta */
  ka::RW<ka::range2d<double> > /* C */
>{};
template<>
struct TaskBodyCPU<TaskDSYRK> {
  void operator()(
    CBLAS_UPLO uplo,
    CBLAS_TRANSPOSE trans,
    double alpha,
    ka::range2d_r <double>  A, 
    double beta,
    ka::range2d_rw<double> C 
  )
  {
    const int n     = A->dim(0); 
    const int k     = A->dim(1); // eq. to Akj.rows();
    const int lda   = A->lda();
    const double* const a = A->ptr();

    const int ldc   = C->lda();
    double* const c = C->ptr();

    cblas_dsyrk
    (
      CblasRowMajor, uplo, trans,
      n, k, alpha, a, lda, beta, c, ldc
    );

  }
};


/* DTRSM
*/
struct TaskDTRSM: public ka::Task<7>::Signature<
      CBLAS_SIDE,                  /* side */
      CBLAS_UPLO,                  /* uplo */
      CBLAS_TRANSPOSE,             /* transA */
      CBLAS_DIAG,                  /* diag */
      double,                      /* alpha */
      ka::R<ka::range2d<double> >, /* A */
      ka::RW<ka::range2d<double> > /* B */
>{};
template<>
struct TaskBodyCPU<TaskDTRSM> {
  void operator()( 
    CBLAS_SIDE             side,
    CBLAS_UPLO             uplo,
    CBLAS_TRANSPOSE        transA,
    CBLAS_DIAG             diag,
    double                 alpha,
    ka::range2d_r <double> A, 
    ka::range2d_rw<double> C
  )
  {
    const double* const a = A->ptr();
    const int lda = A->lda();

    double* const c = C->ptr();
    const int ldc = C->lda();

    const int n = C->dim(0);
    const int k = (transA == CblasNoTrans ? A->dim(1) : A->dim(0) );

    cblas_dtrsm
    (
      CblasRowMajor, side, uplo, transA, diag,
      n, k, alpha, a, lda, c, ldc
    );
  }
};


/* =================== CLAPACK routines =================== */

/* Compute inplace LU factorization of A.
*/
struct TaskDGETRF: public ka::Task<2>::Signature<
  ka::RW<ka::range2d<double> >, /* A */
  ka::W<ka::range1d <int> >  /* pivot */
>{};
template<>
struct TaskBodyCPU<TaskDGETRF> {
  void operator()( 
    ka::range2d_rw<double> A, 
    ka::range1d_w<int> piv
  )
  {
    const int m        = A->dim(0); 
    const int n        = A->dim(0); 
    const int lda      = A->lda();
    double* const a    = A->ptr();
    int* const ipiv    = piv->ptr();

    clapack_dgetrf(
      CblasRowMajor, m, n, a, lda, ipiv
    );
  }
};


/* Compute inplace LLt factorization of A, ie L such that A = L * Lt
   with L lower triangular.
*/
struct TaskDPOTRF: public ka::Task<2>::Signature<
  CBLAS_UPLO,                  /* upper / lower */
  ka::RW<ka::range2d<double> > /* A */
>{};
template<>
struct TaskBodyCPU<TaskDPOTRF> {
  void operator()( 
    CBLAS_UPLO uplo, ka::range2d_rw<double> A 
  )
  {
    const int n     = A->dim(0); 
    const int lda   = A->lda();
    double* const a = A->ptr();
    clapack_dpotrf(
      CblasRowMajor, uplo, n, a, lda
    );
  }
};


/* =================== Misc routines =================== */

/* Compute the norm_infty of A-B 
*/
/* Task Error
 * Compute the || ||2 matrix norm of A-B
 */
struct TaskNorm2 : public ka::Task<3>::Signature<
  ka::W<double>,               /* norm */
  ka::R<ka::range2d<double> >, /* A */
  ka::R<ka::range2d<double> >  /* B */
> {};

template<>
struct TaskBodyCPU<TaskNorm2> {
  void operator() ( 
    ka::pointer_w<double> norm,
    ka::range2d_r<double> A, 
    ka::range2d_r<double> B  
  )
  {
    size_t d0 = A->dim(0);
    size_t d1 = A->dim(1);
    double error = 0.0;
    for (size_t i=0; i < d0; ++i)
    {
      for (size_t j=0; j < d1; ++j)
      {
        double diff = fabs(A(i,j)-B(i,j));
        error += diff*diff;
      }
    }
    *norm = sqrt(error);
  }
};

/********************************************************
 * end: from matrix.h in examples of xkaapi repository
 *******************************************************/



/** Parallel bloc matrix product (lower level)
*/
struct TaskMatProductL2: public ka::Task<3>::Signature<
      ka::R<ka::range2d<double> >, /* A */
      ka::R<ka::range2d<double> >,  /* B */
      ka::RPWP<ka::range2d<double> >   /* C */
>{};

/** Parallel bloc matrix product (upper level)
*/
struct TaskMatProductL1: public ka::Task<3>::Signature<
      ka::R<ka::range2d<double> >, /* A */
      ka::R<ka::range2d<double> >,  /* B */
      ka::RPWP<ka::range2d<double> >   /* C */
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

// multiplies A*B^T and stores it in *this
void multKAAPI(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose);

#endif
#endif
