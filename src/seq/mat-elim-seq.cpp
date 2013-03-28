/**
 * \file   mat-elim-seq.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for sequential Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim-seq.h"

void elimSEQ(Matrix& A, int blocksize) {
  blockElimSEQ(A, 
}

void elimNaiveSEQ(Matrix& A, int blocksize) {
  int m = A.nRows();
  int n = A.nCols(); 
  for (int i = 0; i < m; 
