/***********************************************************************************/ /**
                                                                                       *
                                                                                       * \file SIMD_DEFINITIONS.hpp
                                                                                       *
                                                                                       * \brief Contains macro
                                                                                       *definitions for the intrinsics
                                                                                       *used for the vectorization
                                                                                       *
                                                                                       * \author Wolfgang Hölzl
                                                                                       *(hoelzlw), hoelzlw AT in.tum.de
                                                                                       *
                                                                                       **************************************************************************************/

#pragma once

#ifndef SIMD_DEFINITIONS_H
#define SIMD_DEFINITIONS_H

#include <cmath>
#include <limits>

/*
 * Check whether the file SIMD_TYPES.hpp has been included.
 * This file (SIMD_DEFINITIONS.hpp) needs some macros to be set properly by that file (SIMD_TYPES.hpp).
 */
#ifndef SIMD_TYPES_H
#error "SIMD_DEFINITIONS included without SIMD_TYPES! Never include this file directly! Include it only via SIMD_TYPES!"
#endif /* defined SIMD_TYPES_H */

#if defined VECTOR_SSE4_FLOAT32
/*
 * Map single precision SSE intrinsics
 */
#define ADDV _mm_add_ps
#define SUBV _mm_sub_ps
#define MULV _mm_mul_ps
#define DIVV _mm_div_ps
#define SQRTV _mm_sqrt_ps

#define LOADU _mm_loadu_ps
#define STOREU _mm_storeu_ps
#define SETV_R _mm_set1_ps
#define SETV_I _mm_set1_epi32
#define ZEROV_R _mm_setzero_ps
#define ZEROV_I _mm_setzero_si128

#define MAXV _mm_max_ps
#define MINV _mm_min_ps

#define CMP_LT _mm_cmplt_ps
#define CMP_LE _mm_cmple_ps
#define CMP_GT _mm_cmpgt_ps
#define CMP_GE _mm_cmpge_ps

#define CMP_EQ_I _mm_cmpeq_epi32
#define CMP_EQ_R _mm_cmpeq_ps

#define ANDV_R _mm_and_ps
#define ORV_R _mm_or_ps
#define XORV_R _mm_xor_ps
#define ORV_I _mm_or_si128
#define NOTV_R not_ps
#define NOTV_I not_si128
#define ANDNOTV_R _mm_andnot_ps

#define BLENDV _mm_blendv_ps
#define BLENDV_I(else_part, if_part, mask) \
  CAST_REAL_TO_INT_V(_mm_blendv_ps(CAST_INT_TO_REAL_V(else_part), CAST_INT_TO_REAL_V(if_part), mask))
#define MOVEMASK _mm_movemask_ps
#define SHIFT_LEFT _mm_slli_epi32

#define CAST_INT_TO_REAL_V _mm_castsi128_ps
#define CAST_REAL_TO_INT_V _mm_castps_si128
#define FABS fabs_ps

/*
 * Compute the absolute value of a vector by forcing the sign bit to be zero
 */
inline __m128 fabs_ps(const __m128 x) {
  static const __m128 sign_mask = CAST_INT_TO_REAL_V(_mm_set1_epi32(1 << 31));
  return _mm_andnot_ps(sign_mask, x);
}

/*
 * Bitwise NOT operation for integers
 */
inline __m128i not_si128(const __m128i x) {
  static const __m128i mask = _mm_set1_epi32(~0);
  return CAST_REAL_TO_INT_V(_mm_xor_ps(CAST_INT_TO_REAL_V(mask), CAST_INT_TO_REAL_V(x)));
}

/*
 * Bitwise NOT operation for reals
 */
inline __m128 not_ps(const __m128 x) {
  static const __m128i mask = _mm_set1_epi32(~0);
  return _mm_xor_ps(CAST_INT_TO_REAL_V(mask), x);
}

/*
 * Check, whether a real_vector contains infinity or NaN
 */
inline bool checkVector(const __m128 x) {
  static const real_vector infinity = SETV_R(std ::numeric_limits<float>::infinity());
  return MOVEMASK(ANDV_R(CMP_EQ_R(x, x), NOTV_R(CMP_EQ_R(x, infinity)))) == VECTOR_FULL_MASK;
}
#elif defined VECTOR_SSE4_FLOAT64
/*
 * Map double precision SSE intrinsics
 */
#error "SSE4 with double precision not implemented at the moment"
#elif defined VECTOR_AVX_FLOAT32
/*
 * Map single precision AVX intrinsics
 */
#define ADDV _mm256_add_ps
#define SUBV _mm256_sub_ps
#define MULV _mm256_mul_ps
#define DIVV _mm256_div_ps
#define SQRTV _mm256_sqrt_ps

#define LOADU _mm256_loadu_ps
#define STOREU _mm256_storeu_ps
#define SETV_R _mm256_set1_ps
#define SETV_I _mm256_set1_epi32
#define ZEROV_R _mm256_setzero_ps
#define ZEROV_I _mm256_setzero_si256

#define MAXV _mm256_max_ps
#define MINV _mm256_min_ps

#define CMP_LT(x, y) _mm256_cmp_ps((x), (y), _CMP_LT_OS)
#define CMP_LE(x, y) _mm256_cmp_ps((x), (y), _CMP_LE_OS)
#define CMP_GT(x, y) _mm256_cmp_ps((x), (y), _CMP_GT_OS)
#define CMP_GE(x, y) _mm256_cmp_ps((x), (y), _CMP_GE_OS)

#define CMP_EQ_R(x, y) _mm256_cmp_ps((x), (y), _CMP_EQ_OS)

/*
 * Define test for equality of integers
 * Replace with
 *
 * #define CMP_EQ_I(x, y) _mm256_cmpeq_epi32((x), (y))
 *
 * when running with AVX2
 */
static inline __m256i CMP_EQ_I(const __m256i a, const __m256i b) {
  __m256i              out = ZEROV_I();
  const integer* const p   = reinterpret_cast<const integer*>(&a);
  const integer* const q   = reinterpret_cast<const integer*>(&b);
  integer* const       r   = reinterpret_cast<integer*>(&out);

  for (int i = 0; i < VECTOR_LENGTH; ++i) {
    r[i] = p[i] == q[i] ? 0xFFFFFFFF : 0;
  }

  return out;
}

#define ANDV_R _mm256_and_ps
#define ORV_R _mm256_or_ps
#define XORV_R _mm256_xor_ps
#define ORV_I(x, y) CAST_REAL_TO_INT_V(_mm256_or_ps(CAST_INT_TO_REAL_V(x), CAST_INT_TO_REAL_V(y)))

#define NOTV_R not_ps
#define NOTV_I not_si256
#define ANDNOTV_R _mm256_andnot_ps

#define BLENDV _mm256_blendv_ps
#define BLENDV_I(else_part, if_part, mask) \
  CAST_REAL_TO_INT_V(_mm256_blendv_ps(CAST_INT_TO_REAL_V(else_part), CAST_INT_TO_REAL_V(if_part), mask))
#define MOVEMASK _mm256_movemask_ps

/*
 * Define left shifting for integers
 * Replace with
 *
 * #define SHIFT_LEFT(x, y) _mm256_slli_epi32((x), (y))
 *
 * when running with AVX2
 */
static inline __m256i SHIFT_LEFT(const __m256i x, const integer y) {
  __m256i              out = ZEROV_I();
  const integer* const p   = reinterpret_cast<const integer*>(&x);
  integer* const       q   = reinterpret_cast<integer*>(&out);

  for (int i = 0; i < VECTOR_LENGTH; ++i) {
    q[i] = p[i] << y;
  }

  return out;
}

#define CAST_INT_TO_REAL_V _mm256_castsi256_ps
#define CAST_REAL_TO_INT_V _mm256_castps_si256
#define FABS fabs_ps

/*
 * Compute the absolute value of a vector by forcing the sign bit to be zero
 */
inline __m256 fabs_ps(const __m256 x) {
  static const __m256 sign_mask = CAST_INT_TO_REAL_V(_mm256_set1_epi32(1 << 31));
  return _mm256_andnot_ps(sign_mask, x);
}

/*
 * Bitwise NOT operation for integers
 */
inline __m256i not_si256(const __m256i x) {
  static const __m256i mask = _mm256_set1_epi32(0xFFFFFFFF);
  return CAST_REAL_TO_INT_V(_mm256_xor_ps(CAST_INT_TO_REAL_V(mask), CAST_INT_TO_REAL_V(x)));
}

/*
 * Bitwise NOT operation for reals
 */
inline __m256 not_ps(const __m256 x) {
  static const __m256i mask = _mm256_set1_epi32(0xFFFFFFFF);
  return _mm256_xor_ps(CAST_INT_TO_REAL_V(mask), x);
}

/*
 * Check, whether a real_vector contains infinity or NaN
 */
inline bool checkVector(const __m256 x) {
  static const real_vector infinity = SETV_R(std ::numeric_limits<float>::infinity());
  return MOVEMASK(ANDV_R(CMP_EQ_R(x, x), NOTV_R(CMP_EQ_R(x, infinity)))) == VECTOR_FULL_MASK;
}
#elif defined VECTOR_AVX_FLOAT64
/*
 * Map double precision AVX intrinsics
 */
#error "AVX with double precision not implemented at the moment"
#else /* no vectorization type defined */
/*
 * No vectorization demanded.
 * Do nothing, but inform the user
 */
#pragma message "SIMD-Definitions included, but no Vector-Type defined."
#endif

#endif /* #ifndef SIMD_DEFINITIONS_H */
