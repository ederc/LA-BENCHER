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
    std::cout << "Generate new random matrix." << std::endl;
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
/*
void genMatrix(int rows=0, int cols=0) {
  uint32 m, n;
  bool cmp;
  if (rows == 0 || cols == 0) {
    std::cout << "Generate new random matrix with entries of type float." << std::endl;
    std::cout << "Number of rows (<2^32): ";
    std::cin >> m;
    std::cout << "Number of cols (<2^32): ";
    std::cin >> n;
    std::cout << "Check if matrix is stored correctly? (1=yes, 0=no)  ";
    std::cin >> cmp;
    genRandom(m, n, cmp, 1);
  }
  std::cout << "Matrix generated." << std::endl;
}
*/
