/**
 * \file   mat-elim-seq.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for sequential Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim-seq.h"

mat negInverseModP(mat a, mat p) {
  mat r1  = a;
  mat r2  = p;
  mat rs  = 0;
  mat u1  = 1;
  mat u2  = 0;
  mat us  = 0;
  mat q   = 0;
  while (r2 != 0) {
    q   = r1 / r2;
    rs  = r1;
    us  = u1;
    r1  = r2;
    u1  = u2;
    r2  = rs - q * r2;
    u2  = us - q * u2;
  }
  return p - u1;
}

void cleanUpModP(Matrix& A, mat p) {
  for (int i = 0; i < A.entries.size(); ++i)
    A.entries[i] %= p;
}

void elimSEQ(Matrix& A, int blocksize) {
  //blockElimSEQ(A, 
}

void elimNaiveSEQModP(Matrix& A, int blocksize) {
  // biggest prime < 2^16
  mat prime     = 65521;
  int m         = A.nRows();
  int n         = A.nCols(); 
  // if m > n then only n eliminations are possible
  int boundary  = (m > n) ? n : m;
  mat inv, mult;
  for (int i = 0; i < boundary; ++i) {
    inv  = negInverseModP(A(i,i), prime);
    //std::cout << "A(" << i << "," << i << ") " << A(i,i) << std::endl;
    //std::cout << "inv  " << inv << std::endl;
    for (int j = i+1; j < m; ++j) {
      //std::cout << "A(" << j << "," << i << ") " << A(j,i) << std::endl;
      mult  = A(j,i) * inv;
      for (int k = i; k < n; ++k) {
        A(j,k) += (A(i,k) * mult % prime);
      }
    }
  }
  cleanUpModP(A, prime);
}
