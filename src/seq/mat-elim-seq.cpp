/**
 * \file   mat-elim-seq.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for sequential Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "../modular.h"
#include "mat-elim-seq.h"

#define F4RT_DBG  0

void cleanUpModP(Matrix& A, uint64 p) {
  for (int i = 0; i < A.entries.size(); ++i)
    A.entries[i] %= p;
}

void elimSEQ(Matrix& A, int blocksize) {
  //blockElimSEQ(A, 
}

void elimNaiveSEQModP(Matrix& A, int blocksize, uint64 prime) {
  int l;
  int m         = A.nRows();
  int n         = A.nCols(); 
  // if m > n then only n eliminations are possible
  int boundary  = (m > n) ? n : m;
  mat inv, mult;
  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Naive Gaussian Elimination" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();
  for (int i = 0; i < boundary; ++i) {
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
        for (int it = i; it < n; ++it)
          A.entries[it+l*n] = A.entries[it+i*n];
        for (int it = i; it < n; ++it) {
          std::cout << "it - " << it << " i - " << i << std::endl;
          std::cout << "A.entries " << A.entries[it+i*n] << std::endl;
          std::cout << "tempRow " << tempRow[it-i] << std::endl;
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
    for (int j = i+1; j < m; ++j) {
#if F4RT_DBG
      std::cout << "A(" << j << "," << i << ") " << A(j,i) << std::endl;
#endif
      mult  = (A(j,i) * inv) % prime;
      for (int k = i; k < n; ++k) {
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
  //cleanUpModP(A, prime);
  //A.print();
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Raw sequential" << std::endl;
  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows();
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
