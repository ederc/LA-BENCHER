/**
 * \file   intrinsics.h
 * \author Christian Eder ( ederc@mathematik.uni-kl.de )
 * \date   February 2013
 * \brief  General header file for intrinsic stuff.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_INTRINSICS_H
#define F4RT_INTRINSICS_H

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

#include "f4rt-config.h"

#if defined(__F4RT_HAVE_AVX)
#include <immintrin.h>
#elif defined(__F4RT_HAVE_AES)
#include <wmmintrin.h>
#elif defined(__F4RT_HAVE_SSE4A)
#include <ammintrin.h>
#elif defined(__F4RT_HAVE_SSE4_2)
#include <nmmintrin.h>
#elif defined(__F4RT_HAVE_SSE4_1)
#include <smmintrin.h>
/*
#elif defined(__F4RT_HAVE_SSSE3)
#include <tmmintrin.h>
*/
#elif defined(__F4RT_HAVE_SSE3)
#include <emmintrin.h>
#elif defined(__F4RT_HAVE_SSE2)
#include <emmintrin.h>
#elif defined(__F4RT_HAVE_SSE)
#include <xmmintrin.h>
#elif defined(__F4RT_HAVE_MMX)
#include <mmintrin.h>
#endif

#define X_F4RT

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __F4RT_HAVE_SSE2 
/**
 * \fn  static inline void _add_mm128_epi32(
 *        uint32_t *dst,
 *        uint32_t const *src
 *      )
 *
 * \brief Adds 32 bit unsigned integers stored in a 128 bit vector using 
 * MMX/SSE intrinsics. 
 *
 * \param dst first vector and destination of the result
 *
 * \param src second vector, added to dst
 *
 * \note Both vectors need to be 16 bit aligned, otherwise correctness cannot be
 * guaranteed.
 */
static inline void _add_mm128_epi32(
    uint32_t * dst,
    uint32_t * const src
) {
  /*
  const unsigned int __attribute__ ((aligned(16))) src1[4] = { 
    src1_1,
    src1_2,
    src1_3,
    src1_4
  };
  unsigned int __attribute__ ((aligned(16))) src2[4] = { 
    src2_1,
    src2_2,
    src2_3,
    src2_4
  };
*/
  //unsigned int __attribute__ ((aligned(16))) dst[4];
  __m128i __src  = _mm_load_si128((__m128i*)src);
  __m128i __dst  = _mm_load_si128((__m128i*)dst);
  __dst   = _mm_add_epi32(__src,__dst);
  _mm_store_si128((__m128i*)dst, __dst);
}

/*
inline void _add( uint16_t * dst, uint16_t const * src, size_t n )
{
  for( uint16_t const * end( dst + n ); dst != end; dst+=8, src+=8 )
  {
    __m128i _s = _mm_load_si128( (__m128i*) src );
    __m128i _d = _mm_load_si128( (__m128i*) dst );

    _d = _mm_add_epi16( _d, _s );

    _mm_store_si128( (__m128i*) dst, _d );
  }
}
*/


#endif

#ifdef __cplusplus
}
#endif
#endif
