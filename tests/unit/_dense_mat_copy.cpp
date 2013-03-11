/**
 * \file   _dense_mat_copy.cpp
 * \author Christian Eder ( ederc@mathematik.uni-kl.de )
 * \date   March 2013
 * \brief  Unit test of matrix copying.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "f4rt-config.h"
#include "matrix.h"
#include <assert.h>
int main() {

  // generate matrix of size 20*30 with trivial entries of type float
  Matrix A(20,30,1);
  
  // copy A to B
  Matrix B;
  B.copy(A);

  assert(check(A,B,1) == 0);

  // clear memory
  A.clear();
  B.clear();

  // generate matrix of size 20*30 with random entries of type float
  Matrix C(20,30);
  C.entries.resize(20*30);
  srand(time(NULL));
  for (size_t i = 0; i < 20*30; ++i) {
    C.entries[i] = getRandom();
  }
  
  // copy C to D
  Matrix D;
  D.copy(C);

  assert(check(C,D,1) == 0);

  // clear memory
  C.clear();
  D.clear();


  return 0;
}
