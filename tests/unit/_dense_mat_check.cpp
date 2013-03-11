/**
 * \file   _dense_mat_check.cpp
 * \author Christian Eder ( ederc@mathematik.uni-kl.de )
 * \date   March 2013
 * \brief  Unit test of matrix check.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "f4rt-config.h"
#include "matrix.h"
#include <assert.h>
int main() {

  // generate matrices of size 20*30 with trivial entries of type float
  Matrix A(20,30,1);
  Matrix B(20,30,2);
  Matrix C(20,30,2);

  // test check function
  assert(check(A,B,1) == 1);
  assert(check(A,C,1) == 1);
  assert(check(B,C,1) == 0);

  // clear memory
  A.clear();
  B.clear();
  C.clear();

  return 0;
}
