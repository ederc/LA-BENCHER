/**
 * \file   _mullo_mm128_epi32.c
 * \author Christian Eder ( ederc@mathematik.uni-kl.de )
 * \date   February 2013
 * \brief  Unit test of _mullo_mm128_epi32.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "f4rt-config.h"
#include "intrinsics.h"
#include <assert.h>
int main() {
  uint32_t __attribute__ ((aligned(16))) v[4] = { 65537,233,12365,34};
  uint32_t __attribute__ ((aligned(16))) w[4] = { 65536,2,2,31};
  _mullo_mm128_epi32(v,w);
  assert(
    v[0] == 65536   &&
    v[1] == 466     &&
    v[2] == 24730   &&
    v[3] == 1054
  );
  return 0;
}
