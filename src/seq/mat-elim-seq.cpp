/**
 * \file   mat-elim-seq.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for sequential Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim-seq.h"

#define F4RT_DBG  0

void elimCoSEQBaseModP( mat *M, const uint32 k1, const uint32 i1, 
                        const uint32 j1, const uint32 rows, const uint32 cols,
                        uint64 size, uint64 prime, mat *neg_inv_piv) {
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
    const uint64 istart  = (k1 == i1) ? k+1 : 0;
    for (uint64 i = istart; i < size; i++) {
      const mat tmp = (M[k+k1+(i1+i)*cols] * inv_piv) % prime;
      // if the pivots are in the same column part of the matrix as Mmdf then we can
      // always start at the next column (k+1), otherwise we need to start at
      // column 0
      const uint64 jstart  = (k1 == j1) ? k+1 : 0;
      for (uint64 j = jstart; j < size; j++) {
  	    M[(j1+j)+(i1+i)*cols]  +=  M[(j1+j)+(k1+k)*cols] * tmp;
  	    M[(j1+j)+(i1+i)*cols]  %=  prime;
      }
    }
  }
}

void elimCoSEQBlockModP(mat *M, const uint32 k1, const uint32 k2, 
                        const uint32 i1, const uint32 i2,
		                    const uint32 j1, const uint32 j2, 
		                    const uint32 rows, const uint32 cols, 
                        uint64 size, uint64 prime, mat *neg_inv_piv) {
  if (i2 <= k1 || j2 <= k1) 
    return;
  // 
  //if (size <= 2) {
  if (size <= __F4RT_CPU_L1_CACHE) {
    elimCoSEQBaseModP (M, k1, i1, j1, rows, cols, size, prime, neg_inv_piv);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    elimCoSEQBlockModP( M, k1, km, i1, im, j1, jm, rows, cols, size,
                        prime, neg_inv_piv);
    elimCoSEQBlockModP( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
                        prime, neg_inv_piv);
    elimCoSEQBlockModP( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
                        prime, neg_inv_piv);
    elimCoSEQBlockModP( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
                        prime, neg_inv_piv);
    elimCoSEQBlockModP( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
                        prime, neg_inv_piv);
    elimCoSEQBlockModP( M, km+1, k2, im+1, i2, j1, jm, rows, cols, size,
                        prime, neg_inv_piv);
    elimCoSEQBlockModP( M, km+1, k2, i1, im, jm+1, j2, rows, cols, size,
                        prime, neg_inv_piv);
    elimCoSEQBlockModP( M, km+1, k2, i1, im, j1, jm, rows, cols, size,
                        prime, neg_inv_piv);
  }
}

void D1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime,
        mat *neg_inv_piv, uint32 blocksize) {
  if (i2 <= k1 || j2 <= k1) 
    return;

  if (size <= blocksize) {
    elimCoSEQBaseModP (M, k1, i1, j1, rows, cols, size, prime, neg_inv_piv);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // parallel - start
    // X11
    D1( M, k1, km, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X12
    D1( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X21
    D1( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X22
    D1( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end
    
    // parallel - start
    // X11
    D1( M, km+1, k2, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X12
    D1( M, km+1, k2, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X21
    D1( M, km+1, k2, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X22
    D1( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end
  }
}

void C1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime,
        mat *neg_inv_piv, uint32 blocksize) {
  if (i2 <= k1 || j2 <= k1) 
    return;

  if (size <= blocksize) {
    elimCoSEQBaseModP (M, k1, i1, j1, rows, cols, size, prime, neg_inv_piv);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // parallel - start
    // X11
    C1( M, k1, km, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X21
    C1( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end

    // parallel - start
    // X12
    D1( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X22
    D1( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end
    
    // parallel - start
    // X12
    C1( M, km+1, k2, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X22
    C1( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end
    
    // parallel - start
    // X11
    D1( M, km+1, k2, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X12
    D1( M, km+1, k2, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end
  }
}

void B1(mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime,
        mat *neg_inv_piv, uint32 blocksize) {
  if (i2 <= k1 || j2 <= k1) 
    return;

  if (size <= blocksize) {
    elimCoSEQBaseModP (M, k1, i1, j1, rows, cols, size, prime, neg_inv_piv);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // parallel - start
    // X11
    B1( M, k1, km, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X12
    B1( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end

    // parallel - start
    // X21
    D1( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X22
    D1( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end
    
    // parallel - start
    // X21
    B1( M, km+1, k2, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X22
    B1( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end
    
    // parallel - start
    // X11
    D1( M, km+1, k2, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // X12
    D1( M, km+1, k2, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end

  }
}

void A( mat *M, const uint32 k1, const uint32 k2, 
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2, 
		    const uint32 rows, const uint32 cols, 
        uint64 size, uint64 prime,
        mat *neg_inv_piv, uint32 blocksize) {
  if (i2 <= k1 || j2 <= k1) 
    return;
  // 
  //if (size <= 2) {
  if (size <= blocksize) {
    elimCoSEQBaseModP (M, k1, i1, j1, rows, cols, size, prime, neg_inv_piv);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // forward step

    A(M, k1, km, i1, im, j1, jm, rows, cols, size,
      prime, neg_inv_piv, blocksize);
    // parallel - start
    B1( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    C1( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    // parallel - end
    D1( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, blocksize);
    
    // backward step
    
    A(M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
      prime, neg_inv_piv, blocksize);
  }
}

void elimCoSEQModPOld(Matrix& A, uint32 blocksize, uint64 prime) {
  uint32 m          = A.nRows();
  uint32 n          = A.nCols();
  // if m > n then only n eliminations are possible
  uint32 boundary   = (m > n) ? n : m;
  mat *a_entries    = A.entries.data();
  mat *neg_inv_piv  =   (mat *)calloc(boundary, sizeof(mat));
  a_entries[0]      %=  prime;
  neg_inv_piv[0]    =   negInverseModP(a_entries[0], prime); 
  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Cache-oblivious Gaussian Elimination without pivoting -- ";
  std::cout << " old version ( not parallelizable )" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();
  
  // computation of blocks
  elimCoSEQBlockModP( a_entries, 0, boundary-1, 0, m-1, 0, n-1, m, n, boundary,
                      prime, neg_inv_piv);
  
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Raw sequential" << std::endl;
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
  std::cout << "# Threads:        " << 1 << std::endl;
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

void elimCoSEQModP(Matrix& M, uint32 blocksize, uint64 prime) {
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
  A(a_entries, 0, boundary-1, 0, m-1, 0, n-1, m, n, boundary,
                              prime, neg_inv_piv, blocksize);
  
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Raw sequential" << std::endl;
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
  std::cout << "# Threads:        " << 1 << std::endl;
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

void elimNaiveSEQModP(Matrix& A, uint32 blocksize, uint64 prime) {
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
    A(i,i) %= prime;
#if F4RT_DBG
    std::cout << "!! A(" << i << "," << i << ") " << A(i,i) << std::endl;
    std::cout << "A(" << i << "," << i << ") " << A(i,i) % prime << std::endl;
#endif
    /*
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
    */
    inv  = negInverseModP(A(i,i), prime);
#if F4RT_DBG
    std::cout << "inv  " << inv << std::endl;
#endif
    for (uint32 j = i+1; j < m; ++j) {
#if F4RT_DBG
      std::cout << "A(" << j << "," << i << ") " << A(j,i) << std::endl;
#endif
      mult  = (A(j,i) * inv) % prime;
      for (uint32 k = i; k < n; ++k) {
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
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Raw sequential" << std::endl;
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
  std::cout << "# Threads:        " << 1 << std::endl;
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

void elimNaiveSEQModPPivot(Matrix& A, uint32 blocksize, uint64 prime) {
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
    for (uint32 j = i+1; j < m; ++j) {
#if F4RT_DBG
      std::cout << "A(" << j << "," << i << ") " << A(j,i) << std::endl;
#endif
      mult  = (A(j,i) * inv) % prime;
      for (uint32 k = i+1; k < n; ++k) {
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
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Raw sequential" << std::endl;
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
  std::cout << "# Threads:        " << 1 << std::endl;
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
