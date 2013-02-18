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

#ifdef __F4RT_HAVE_AVX
#include <immintrin.h>
#endif
#ifdef __F4RT_HAVE_AES
#include <wmmintrin.h>
#endif
#ifdef __F4RT_HAVE_SSE4A
#include <ammintrin.h>
#endif
#ifdef __F4RT_HAVE_SSE42
#include <nmmintrin.h>
#endif
#ifdef __F4RT_HAVE_SSE41
#include <smmintrin.h>
#endif
#ifdef __F4RT_HAVE_SSSE3
#include <tmmintrin.h>
#endif
#ifdef __F4RT_HAVE_SSE3
#include <emmintrin.h>
#endif
#ifdef __F4RT_HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef __F4RT_HAVE_SSE
#include <xmmintrin.h>
#endif
#ifdef __F4RT_HAVE_MMX
#include <mmintrin.h>
#endif

#define X_F4RT

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __F4RT_HAVE_SSE41
/**
 * \fn  static inline void _add_mm128_epi32(
 *        uint32_t *dst,
 *        uint32_t const *src
 *      )
 *
 * \brief Adds 16 bit unsigned integers stored in a 128 bit vector using 
 * MMX/SSE intrinsics. Those integers are stored in vectors of uint32_t, for the
 * sake of forbidding overflows.
 *
 * \param dst first vector and destination of the result
 *
 * \param src second vector, added to dst
 *
 * \note Both vectors need to be 16 bit aligned, otherwise correctness cannot be
 * guaranteed.
 * Moreover, the integers in the vectors need to be <=2^31, otherwise an
 * overflow could occur: 2^31 + 2^31 = 2 * 2^31 = 2^32.
 */
static inline void _add_mm128_epi32(
    uint32_t * dst,
    uint32_t * const src
) {
  __m128i __src  = _mm_load_si128((__m128i*)src);
  __m128i __dst  = _mm_load_si128((__m128i*)dst);
  __dst   = _mm_add_epi32(__src,__dst);
  _mm_store_si128((__m128i*)dst, __dst);
}

/**
 * \fn  static inline void _mul_mm128_epi32(
 *        int32_t *dst,
 *        int32_t const *src
 *      )
 *
 * \brief Multiplies 32 bit signed integers stored in a 128 bit vector using 
 * MMX/SSE intrinsics. Those integers are stored in vectors of int32_t, for the
 * sake of forbidding overflows.
 *
 * \param dst first vector and destination of the result
 *
 * \param src second vector, multiplied to dst
 *
 * \note Both vectors need to be 16 bit aligned, otherwise correctness cannot be
 * guaranteed.
 * Moreover, the integers in the vectors need to be <=2^16, otherwise an
 * overflow could occur: 2^16 * 2^16 = 2^(16+16) = 2^32.
 */
static inline void _mul_mm128_epi32(
    int32_t * dst,
    int32_t * const src
) {
  __m128i __src  = _mm_load_si128((__m128i*)src);
  __m128i __dst  = _mm_load_si128((__m128i*)dst);
  __dst   = _mm_mul_epi32(__src,__dst);
  _mm_store_si128((__m128i*)dst, __dst);
}

/**
 * \fn  static inline void _mullo_mm128_epi32(
 *        uint32_t *dst,
 *        uint32_t const *src
 *      )
 *
 * \brief Multiplies 32 bit unsigned integers stored in a 128 bit vector using 
 * MMX/SSE intrinsics. Those integers are stored in vectors of uint32_t, for the
 * sake of forbidding overflows.
 *
 * \param dst first vector and destination of the result
 *
 * \param src second vector, multiplied to dst
 *
 * \note Both vectors need to be 16 bit aligned, otherwise correctness cannot be
 * guaranteed.
 * Moreover, the integers in the vectors need to be <=2^16, otherwise an
 * overflow could occur: 2^16 * 2^16 = 2^(16+16) = 2^32.
 */
static inline void _mullo_mm128_epi32(
    uint32_t * dst,
    uint32_t * const src
) {
  __m128i __src  = _mm_load_si128((__m128i*)src);
  __m128i __dst  = _mm_load_si128((__m128i*)dst);
  __dst   = _mm_mullo_epi32(__src,__dst);
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
