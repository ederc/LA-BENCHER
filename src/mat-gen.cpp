/**
 * \file   mat-gen.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   February 2013
 * \brief  Source file for matrix generation.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-gen.h"

void genMatrix(int rows=0, int cols=0) {
  uint32 m, n;
  bool cmp;
  Matrix A;
  if (rows == 0 || cols == 0) {
    std::cout << "Generate new random matrix with entries of type float." << std::endl;
    std::cout << "Number of rows (<2^32): ";
    std::cin >> m;
    std::cout << "Number of cols (<2^32): ";
    std::cin >> n;
    std::cout << "Check if matrix is stored correctly? (1=yes, 0=no)  ";
    std::cin >> cmp;
    A.generateRandomMatrix(m,n,cmp,1);
  } else {
    A.generateRandomMatrix(rows,cols,0,0);
  }
  A.clear();
  std::cout << "Matrix generated." << std::endl;
}

void genMatrixNew(int rows=0, int cols=0) {
  uint32 m, n, t;
  bool cmp;
  if (rows == 0 || cols == 0) {
    std::cout << "Generate new random matrix with entries of type float." << std::endl;
    std::cout << "Number of rows (<2^32): ";
    std::cin >> m;
    std::cout << "Number of cols (<2^32): ";
    std::cin >> n;
    std::cout << "Type of entries (uint16 = 1, uint32 = 2, float = 3): ";
    std::cin >> t;
    std::cout << "Check if matrix is stored correctly? (1=yes, 0=no)  ";
    std::cin >> cmp;
    if (t == 1) {
      matrix<uint16> M;
      genRandom(M, m, n, t, cmp, 1);
    }
    if (t == 2) {
      matrix<uint32> M;
      genRandom(M, m, n, t, cmp, 1);
    }
    if (t == 3) {
      matrix<float> M;
      genRandom(M, m, n, t, cmp, 1);
    }
  }
  std::cout << "Matrix generated." << std::endl;
}
