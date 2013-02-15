/**
 * \file   intrinsics.c
 * \author Christian Eder ( ederc@mathematik.uni-kl.de )
 * \date   February 2013
 * \brief  General source file for intrinsic stuff.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <intrinsics.h>

int main() {
  unsigned int __attribute__ ((aligned(16))) v[4] = { 234,233,12365,34234};
  unsigned int __attribute__ ((aligned(16))) w[4] = { 4658,1,0,312121};
  unsigned int __attribute__ ((aligned(16))) x[4];
  _add_mm128_epi32(v,w);
  printf("%d - %d - %d - %d\n", v[0],v[1],v[2],v[3]);
  return 0;
}
