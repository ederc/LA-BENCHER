/**
 * \file   mat-elim-tools.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for GEP tools.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim-tools.h"

double countGEPFlops(uint32 m, uint32 n) {
  uint32 boundary = m > n ? n : m;
  double res = 0;
  for (uint32 i = 1; i <= boundary; ++i) {
    res +=  (2*(n-boundary)+1)*(m-boundary-1);
  }
  std::cout << "flops - " << std::setprecision(4) << res << std::endl;
  return res;
}

void cleanUpModP(Matrix& A, uint64 p) {
  for (uint64 i = 0; i < A.entries.size(); ++i)
    A.entries[i] %= p;
}

mat negInverseModP(mat a, uint64 p) {
  // we do two turns of the extended Euclidian algorithm per
  // loop. Usually the sign of x changes each time through the loop,
  // but we avoid that by representing every other x as its negative,
  // which is the value minusLastX. This way no negative values show
  // up.
  mat b           = p;
  mat minusLastX  = 0;
  mat x           = 1;
  while (true) {
    // 1st turn
    if (a == 1)
      break;
    const mat firstQuot =   b / a;
    b                   -=  firstQuot * a;
    minusLastX          +=  firstQuot * x;

    // 2nd turn
    if (b == 1) {
      x = p - minusLastX;
      break;
    }
    const mat secondQuot  =   a / b;
    a                     -=  secondQuot * b;
    x                     +=  secondQuot * minusLastX;
  }
  return p - x;
}
