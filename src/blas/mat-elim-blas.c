/**
 * \file   mat-elim-blas.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for Gaussian Elimination with BLAS.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim-blas.h"
#ifdef __F4RT_HAVE_PLASMA
#include "plasma.h"
#include "core_blas.h"
#include "quark.h"
#endif
#ifdef __F4RT_HAVE_LAPACKE
#include "lapacke.h"
#endif
#ifdef __F4RT_HAVE_CBLAS
#include "cblas.h"
#endif

#define F4RT_DBG  0

#if defined(__F4RT_HAVE_LAPACK) && defined(__F4RT_HAVE_BLAS)

#ifdef __cplusplus
extern "C" {
#endif
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
#ifdef __cplusplus
}
#endif
void elimBLAS(double *M, uint32 rows, uint32 cols, int nthrds, uint32 blocksize, uint64 prime) {
/*
  uint32 m    = rows;
  uint32 n    = cols;
  uint64 dim  = m*n;
  // if m > n then only n eliminations are possible
  //uint32 boundary  = (m > n) ? n : m;
  //double *ipiv  = (double *)malloc(dim * sizeof(double));
  int errorHandler;
  struct timeval start, stop;
  clock_t cStart, cStop;
  printf("Tiled Gaussian Elimination without pivoting\n");
  gettimeofday(&start, NULL);
  cStart  = clock();

  // ATLAS implementation - start
  //errorHandler = ATL_dgetrf(CblasRowMajor, m, n, M, m, ipiv);
  // ATLAS implementation - end

  // LAPACK implementation - start
  //dgetrf_(&m, &n, M, &m, ipiv, &errorHandler);
  // LAPACK implementation - end
  
  // PLASMA implementation - start
  //int *ipiv;
  //PLASMA_desc *L;
  //PLASMA_Init(nthrds);
  //PLASMA_Init_Affinity(nthrds, NULL);
  //errorHandler = PLASMA_Alloc_Workspace_dgetrf_incpiv(m, n, &L, &ipiv);
  //errorHandler = PLASMA_dgetrf_incpiv(m, n, M, m, L, ipiv);
  //PLASMA_Finalize();
  //printf("INFO %d\n", errorHandler);
  // PLASMA implementation - end

  gettimeofday(&stop, NULL);
  cStop = clock();
  printf("---------------------------------------------------\n");
  printf("Method:           OpenBLAS\n");
  // compute FLOPS:
  double flops = countGEPFlopsNoPrime(m, n, prime);
  float epsilon = 0.0000000001;
  double realtime = ((stop.tv_sec - start.tv_sec) * 1e6 +
                    (stop.tv_usec - start.tv_usec)) / 1e6;
  double cputime  = (double)((cStop - cStart)) / CLOCKS_PER_SEC;
  // get digits before decimal point of cputime (the longest number) and setw
  // with it: digits + 1 (point) + 4 (precision)
  double ratio = cputime/realtime;
  printf("# Threads:        %d\n", nthrds);
  printf("Block size:       %u\n", blocksize);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("Real time:        %.4f sec\n", realtime);
  printf("CPU time:         %.4f sec\n", cputime);
  if (cputime > epsilon)
    printf("CPU/real time:    %.4f\n", ratio);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("GFLOPS/sec:       %.4f\n", flops / (1000000000 * realtime));
  printf("---------------------------------------------------\n");
*/
}
#endif
