/**
 * \file   _dense_mult_transpose.cpp
 * \author Christian Eder ( ederc@mathematik.uni-kl.de )
 * \date   March 2013
 * \brief  Unit test of matrix transpose.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "f4rt-config.h"
#include "matrix.h"
#include <assert.h>
int main() {

  // generate a matrix of size 20*30 with random entries of type float
  Matrix A(20,30);
  A.entries.resize(20*30);
  srand(time(NULL));
  for (size_t i = 0; i < 20*30; ++i) {
    A.entries[i] = getRandom();
  }

  // keep entries vector
  std::vector<float> tmp = A.entries;

  // transpose A
  A.transpose();

  assert(A.nRows() == 30);
  assert(A.nCols() == 20);

  for (size_t i=0; i<30; ++i) {
    for (size_t j=0; j<20; ++j) {
      assert(A.entries[j+i*20] = tmp[i+j*30]);
    }
  }
  A.clear();
  tmp.clear();
  return 0;
}
