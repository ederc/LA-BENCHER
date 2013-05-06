/**
 * \file   mat-elim-tools.c
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for GEP tools.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim-tools.h"

uint32 bitlog(uint32 x) {
  uint8    b;
  uint32 res;

  if (x <=  8 ) { /* Shorten computation for small numbers */
    res = 2 * x;
  } else {
    b = 15; /* Find the highest non zero bit in the input argument */
    while ((b > 2) && ((int32)x > 0)) {
      --b;
      x <<= 1;
    }
    x &= 0x7000;
    x >>= 12;

    res = x + 8 * (b - 1);
  }

  return res;
}


double countGEPFlops(uint32 m, uint32 n, uint64 prime) {
  uint32 boundary = m > n ? n : m;
  double logp = (double) bitlog(prime);
  double res = 0;
  for (uint32 i = 1; i <= boundary; ++i) {
    res +=  (2*(n-i)+1+logp)*(m-i) + logp * (n-i) * (m-i);
  }
  return res;
}


double countGEPFlopsNoPrime(uint32 m, uint32 n, uint64 prime) {
  uint32 boundary = m > n ? n : m;
  double res = 0;
  for (uint32 i = 1; i <= boundary; ++i) {
    res +=  (2*(n-i)+1)*(m-i);
  }
  return res;
}

/*
void cleanUpModP(Matrix& A, uint64 p) {
  for (uint64 i = 0; i < A.entries.size(); ++i)
    A.entries[i] %= p;
}
*/

mat negInverseModP(mat a, uint64 p) {
  // we do two turns of the extended Euclidian algorithm per
  // loop. Usually the sign of x changes each time through the loop,
  // but we avoid that by representing every other x as its negative,
  // which is the value minusLastX. This way no negative values show
  // up.
  mat b           = p;
  mat minusLastX  = 0;
  mat x           = 1;
  while (1) {
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
