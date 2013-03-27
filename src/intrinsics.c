/**
 * \file   intrinsics.c
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   February 2013
 * \brief  General source file for intrinsic stuff.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include <intrinsics.h>

int main() {
  uint32_t __attribute__ ((aligned(16))) v[4] = { 65537,233,12365,34};
  uint32_t __attribute__ ((aligned(16))) w[4] = { 65536,2,2,31};
  _mullo_mm128_epi32(v,w);
  printf("%d - %d - %d - %d\n", v[0],v[1],v[2],v[3]);
  return 0;
}
