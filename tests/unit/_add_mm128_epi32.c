/**
 * \file   _add_mm128_epi32.c
 * \author Christian Eder ( ederc@mathematik.uni-kl.de )
 * \date   February 2013
 * \brief  Unit test of _add_mm128_epi32.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "f4rt-config.h"
#include "intrinsics.h"
#include <assert.h>
int main() {
  unsigned int __attribute__ ((aligned(16))) v[4] = { 234,233,12365,34234};
  unsigned int __attribute__ ((aligned(16))) w[4] = { 4658,1,0,312121};
  _add_mm128_epi32(v,w);
  assert(
    v[0] == 4892    &&
    v[1] == 234     &&
    v[2] == 12365   &&
    v[3] == 346355
  );
  return 0;
}
