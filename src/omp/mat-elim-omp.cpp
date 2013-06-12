/**
 * \file   mat-elim-omp.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for Gaussian Elimination using OMP.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim-omp.h"

#define F4RT_DBG  0

#ifdef __F4RT_HAVE_OPENMP
void elimOMP(Matrix& A, uint32 blocksize) {
  //blockElimSEQ(A, 
}

void elimNaiveOMPModP1dOuter(Matrix& A, int nthrds, uint32 blocksize, uint64 prime) {
  int thrdCounter = nthrds;
  uint32 m         = A.nRows();
  uint32 n         = A.nCols(); 
  // if m > n then only n eliminations are possible
  uint32 boundary  = (m > n) ? n : m;
  mat inv, mult;
  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Naive Gaussian Elimination without pivoting" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();
  for (uint32 i = 0; i < boundary; ++i) {
    //A(i,i) %= prime;
#if F4RT_DBG
    std::cout << "!! A(" << i << "," << i << ") " << A(i,i) << std::endl;
    std::cout << "A(" << i << "," << i << ") " << A(i,i) % prime << std::endl;
#endif
    inv  = negInverseModP(A(i,i), prime);
#if F4RT_DBG
    std::cout << "inv  " << inv << std::endl;
#endif
  if (nthrds > 0)
    omp_set_num_threads(nthrds);
  uint32 j, k;
#pragma omp parallel shared(A) private(mult)
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
#pragma omp for schedule(guided,blocksize) private(j, k) collapse(1) 
    for (j = i+1; j < m; ++j) {
#if F4RT_DBG
      std::cout << "A(" << j << "," << i << ") " << A(j,i) << std::endl;
#endif
      mult  = (A(j,i) * inv);// % prime;
      for (k = i+1; k < n; ++k) {
#if F4RT_DBG
        std::cout << "A * mult " << A(i,k)*mult << " - " << (A(i,k)*mult) % prime << " - "
          << (A(i,k)%prime) * (mult % prime) << std::endl;
#endif
        A(j,k) += A(i,k) * mult;
        //A(j,k) %= prime;
#if F4RT_DBG
        std::cout << "A(" << j << "," << k << ") " << A(j,k) << " - " << A(j,k) % prime << std::endl;
#endif
      }
    }
  }
}
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Open MP collapse(1) outer loop" << std::endl;
  // compute FLOPS:
  double flops = countGEPFlops(m, n, prime);
  float epsilon = 0.0000000001;
  double realtime = ((stop.tv_sec - start.tv_sec) * 1e6 + 
                    (stop.tv_usec - start.tv_usec)) / 1e6;
  double cputime  = (double)((cStop - cStart)) / CLOCKS_PER_SEC;
  char buffer[50];
  // get digits before decimal point of cputime (the longest number) and setw
  // with it: digits + 1 (point) + 4 (precision) 
  int digits = sprintf(buffer,"%.0f",cputime);
  double ratio = cputime/realtime;
  std::cout << "# Threads:        " << thrdCounter << std::endl;
  std::cout << "Block size:       " << blocksize << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "Real time:        " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << realtime << " sec" 
    << std::endl;
  std::cout << "CPU time:         " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << cputime
    << " sec" << std::endl;
  if (cputime > epsilon)
    std::cout << "CPU/real time:    " << std::setw(digits+1+4) 
      << std::setprecision(4) << std::fixed << ratio << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "GFLOPS/sec:       " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << flops / (1000000000 * realtime) 
    << std:: endl;
  std::cout << "---------------------------------------------------" << std::endl;
}

void elimNaiveOMPModP1dOuterPivot(Matrix& A, int nthrds, uint32 blocksize,
                                  uint64 prime) {
  int thrdCounter = nthrds;
  uint32 l;
  uint32 m         = A.nRows();
  uint32 n         = A.nCols(); 
  // if m > n then only n eliminations are possible
  uint32 boundary  = (m > n) ? n : m;
  mat inv, mult;
  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Naive Gaussian Elimination with pivoting" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();
  for (uint32 i = 0; i < boundary; ++i) {
    A(i,i) %= prime;
#if F4RT_DBG
    std::cout << "!! A(" << i << "," << i << ") " << A(i,i) << std::endl;
    std::cout << "A(" << i << "," << i << ") " << A(i,i) % prime << std::endl;
#endif
    if (A(i,i) == 0) {
      l = i+1;
#if F4RT_DBG
      std::cout << "l1 " << l << std::endl; 
#endif
      while (l < m && A(l,i) % prime == 0) {
        l++;
#if F4RT_DBG
        std::cout << "l2 " << l << std::endl; 
#endif
      }
      if (l == m) {
        continue;
      } else {
        // swapping
#if F4RT_DBG
        std::cout << "before swapping i = " << i << " with l = " << l << std::endl;
        for (int kk=i; kk < n; ++kk) {
          std::cout << "A(i,"<< kk << ") = " << A(i,kk) << std::endl;
          std::cout << "A(l," << kk << ") = " << A(l,kk) << std::endl;
        }
        std::cout << "n  " << n << std::endl;
        std::cout << "l  " << l << std::endl;
        std::cout << "i  " << i << std::endl;
#endif
        std::vector<mat> tempRow(
            A.entries.begin()+i+(l*n),
            A.entries.begin()+n-1+(l*n)+1);
#if F4RT_DBG
        std::cout << "(n-1)-i " << n-1-i << std::endl;
        std::cout << "sizeof mat " << sizeof(mat) << std::endl;
        std::cout << "sizeof temp " << tempRow.size() << std::endl;
#endif
        for (uint32 it = i; it < n; ++it)
          A.entries[it+l*n] = A.entries[it+i*n];
        for (uint32 it = i; it < n; ++it) {
          A.entries[it+i*n] = tempRow[it-i];
        }
#if F4RT_DBG
        std::cout << "after swapping" << std::endl;
        for (int kk=i; kk < n; ++kk) {
          std::cout << "A(i,"<< kk << ") = " << A(i,kk) << std::endl;
          std::cout << "A(l," << kk << ") = " << A(l,kk) << std::endl;
        }
#endif
        tempRow.clear();
      }
    }
    inv  = negInverseModP(A(i,i), prime);
#if F4RT_DBG
    std::cout << "inv  " << inv << std::endl;
#endif
  if (nthrds > 0)
    omp_set_num_threads(nthrds);
  uint32 j, k;
#pragma omp parallel shared(A) private(mult)
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
#pragma omp for schedule(guided,blocksize) private(j, k) collapse(1) 
    for (j = i+1; j < m; ++j) {
#if F4RT_DBG
      std::cout << "A(" << j << "," << i << ") " << A(j,i) << std::endl;
#endif
      mult  = (A(j,i) * inv) % prime;
      for (k = i+1; k < n; ++k) {
#if F4RT_DBG
        std::cout << "A * mult " << A(i,k)*mult << " - " << (A(i,k)*mult) % prime << " - "
          << (A(i,k)%prime) * (mult % prime) << std::endl;
#endif
        A(j,k) += A(i,k) * mult;
        A(j,k) %= prime;
#if F4RT_DBG
        std::cout << "A(" << j << "," << k << ") " << A(j,k) << " - " << A(j,k) % prime << std::endl;
#endif
      }
    }
  }
}
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Open MP collapse(1) outer loop" << std::endl;
  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = countGEPFlops(m, n, prime);
  float epsilon = 0.0000000001;
  double realtime = ((stop.tv_sec - start.tv_sec) * 1e6 + 
                    (stop.tv_usec - start.tv_usec)) / 1e6;
  double cputime  = (double)((cStop - cStart)) / CLOCKS_PER_SEC;
  char buffer[50];
  // get digits before decimal point of cputime (the longest number) and setw
  // with it: digits + 1 (point) + 4 (precision) 
  int digits = sprintf(buffer,"%.0f",cputime);
  double ratio = cputime/realtime;
  std::cout << "# Threads:        " << thrdCounter << std::endl;
  std::cout << "Block size:       " << blocksize << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "Real time:        " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << realtime << " sec" 
    << std::endl;
  std::cout << "CPU time:         " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << cputime
    << " sec" << std::endl;
  if (cputime > epsilon)
    std::cout << "CPU/real time:    " << std::setw(digits+1+4) 
      << std::setprecision(4) << std::fixed << ratio << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "GFLOPS/sec:       " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << flops / (1000000000 * realtime) 
    << std:: endl;
  std::cout << "---------------------------------------------------" << std::endl;
}


// cache-oblivious stuff

void elimCoOMPBaseModP( mat *M, const uint32 k1, const uint32 i1,
                        const uint32 j1, const uint32 rows, const uint32 cols,
                        uint64 size, uint64 prime, mat *neg_inv_piv, int nthrds) {
  uint64 k;

  for (k = 0; k < size; k++) {
    M[(k1+k)+(k1+k)*cols] %= prime;
    // possibly the negative inverses of the pivots at place (k,k) were already
    // computed in another call. otherwise we need to compute and store it

    if (!neg_inv_piv[k+k1]) {
      if (M[(k1+k)+(k1+k)*cols] != 0) {
        neg_inv_piv[k+k1] = negInverseModP(M[(k1+k)+(k1+k)*cols], prime);
      }
    }
    const mat inv_piv   = neg_inv_piv[k+k1];
    // if the pivots are in the same row part of the matrix as Mmdf then we can
    // always start at the next row (k+1), otherwise we need to start at
    // row 0
    //uint32 i, j;
    const uint64 istart  = (k1 == i1) ? k+1 : 0;

// NOTE: Parallelizing this part of the GEP slows down the computations on my
// machine for 4096 x 4096 by roughly 5%

//#pragma omp parallel
//{
//#pragma omp for schedule(static)
    for (uint64 i = istart; i < size; i++) {
      const mat tmp = (M[k+k1+(i1+i)*cols] * inv_piv);
      //const mat tmp = (M[k+k1+(i1+i)*cols] * inv_piv) % prime;
      // if the pivots are in the same column part of the matrix as Mmdf then we can
      // always start at the next column (k+1), otherwise we need to start at
      // column 0
      const uint64 jstart  = (k1 == j1) ? k+1 : 0;
      for (uint64 j = jstart; j < size; j++) {
        M[(j1+j)+(i1+i)*cols]  +=  M[(j1+j)+(k1+k)*cols] * tmp;
        //M[(j1+j)+(i1+i)*cols]  %=  prime;
      }
    }
//}
  }
}

void D1OMP( mat *M, const uint32 k1, const uint32 k2,
            const uint32 i1, const uint32 i2,
		        const uint32 j1, const uint32 j2,
		        const uint32 rows, const uint32 cols,
            uint64 size, uint64 prime, mat *neg_inv_piv,
            int nthrds, uint32 blocksize) {
  //omp_set_nested(1);
  if (i2 <= k1 || j2 <= k1)
    return;

  if (size <= blocksize) {
    elimCoOMPBaseModP(M, k1, i1, j1, rows, cols, size, prime,
                      neg_inv_piv, nthrds);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X11
    D1OMP( M, k1, km, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X12
    D1OMP( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X21
    D1OMP( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X22
    D1OMP( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X11
    D1OMP( M, km+1, k2, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X12
    D1OMP( M, km+1, k2, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X21
    D1OMP( M, km+1, k2, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X22
    D1OMP( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end
  }
}

void C1OMP(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize) {
  //omp_set_nested(1);
  if (i2 <= k1 || j2 <= k1)
    return;

  if (size <= blocksize) {
    elimCoOMPBaseModP(M, k1, i1, j1, rows, cols, size, prime,
                      neg_inv_piv, nthrds);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X11
    C1OMP( M, k1, km, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X21
    C1OMP( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X12
    D1OMP( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X22
    D1OMP( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X12
    C1OMP( M, km+1, k2, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X22
    C1OMP( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X11
    D1OMP( M, km+1, k2, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X12
    D1OMP( M, km+1, k2, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end
  }
}

void B1OMP(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize) {
  //omp_set_nested(1);
  if (i2 <= k1 || j2 <= k1)
    return;

  if (size <= blocksize) {
    elimCoOMPBaseModP(M, k1, i1, j1, rows, cols, size, prime,
                      neg_inv_piv, nthrds);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X11
    B1OMP( M, k1, km, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X12
    B1OMP( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X21
    D1OMP( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X22
    D1OMP( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X21
    B1OMP( M, km+1, k2, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X22
    B1OMP( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end

    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    // X11
    D1OMP( M, km+1, k2, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    // X12
    D1OMP( M, km+1, k2, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end

  }
}

void AOMP( mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize) {
  //omp_set_nested(1);
  if (i2 <= k1 || j2 <= k1)
    return;

  //if (size <= 2) {
  if (size <= blocksize) {
    elimCoOMPBaseModP(M, k1, i1, j1, rows, cols, size, prime,
                      neg_inv_piv, nthrds);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // forward step

    AOMP(M, k1, km, i1, im, j1, jm, rows, cols, size,
      prime, neg_inv_piv, nthrds, blocksize);
    // parallel - start
#pragma omp parallel sections shared(im,jm,km,prime,M, neg_inv_piv)
{
#pragma omp section
    B1OMP( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
#pragma omp section
    C1OMP( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
}
    // parallel - end
    D1OMP( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);

    // backward step

    AOMP(M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
      prime, neg_inv_piv, nthrds, blocksize);
  }
}

void elimCoOMPModP(Matrix& M, int nthrds, uint32 blocksize, uint64 prime) {
  int thrdCounter   = nthrds;
  omp_set_nested(1);
  if (nthrds > 0)
    omp_set_num_threads(nthrds);
#pragma omp parallel
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
}
  uint32 m          = M.nRows();
  uint32 n          = M.nCols();
  // if m > n then only n eliminations are possible
  uint32 boundary   = (m > n) ? n : m;
  mat *a_entries    = M.entries.data();
  mat *neg_inv_piv  =   (mat *)calloc(boundary, sizeof(mat));
  a_entries[0]      %=  prime;
  neg_inv_piv[0]    =   negInverseModP(a_entries[0], prime);
  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Cache-oblivious Gaussian Elimination without pivoting" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();

  // computation of blocks
  AOMP( a_entries, 0, boundary-1, 0, m-1, 0, n-1, m, n,
        boundary, prime, neg_inv_piv, nthrds, blocksize);

  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Open MP parallel sections" << std::endl;
  // compute FLOPS:
  double flops = countGEPFlops(m, n, prime);
  float epsilon = 0.0000000001;
  double realtime = ((stop.tv_sec - start.tv_sec) * 1e6 +
                    (stop.tv_usec - start.tv_usec)) / 1e6;
  double cputime  = (double)((cStop - cStart)) / CLOCKS_PER_SEC;
  char buffer[50];
  // get digits before decimal point of cputime (the longest number) and setw
  // with it: digits + 1 (point) + 4 (precision)
  int digits = sprintf(buffer,"%.0f",cputime);
  double ratio = cputime/realtime;
  std::cout << "# Threads:        " << thrdCounter << std::endl;
  std::cout << "Block size:       " << blocksize << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "Real time:        " << std::setw(digits+1+4)
    << std::setprecision(4) << std::fixed << realtime << " sec"
    << std::endl;
  std::cout << "CPU time:         " << std::setw(digits+1+4)
    << std::setprecision(4) << std::fixed << cputime
    << " sec" << std::endl;
  if (cputime > epsilon)
    std::cout << "CPU/real time:    " << std::setw(digits+1+4)
      << std::setprecision(4) << std::fixed << ratio << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "GFLOPS/sec:       " << std::setw(digits+1+4)
    << std::setprecision(4) << std::fixed << flops / (1000000000 * realtime)
    << std:: endl;
  std::cout << "---------------------------------------------------" << std::endl;
}
#endif
