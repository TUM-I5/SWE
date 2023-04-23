/***********************************************************************************/ /**
                                                                                       *
                                                                                       * \file SIMD_TYPES.hpp
                                                                                       *
                                                                                       * \brief Defines the length of
                                                                                       *the vectors and the
                                                                                       *corresponding functions
                                                                                       *
                                                                                       * \author Wolfgang Hölzl
                                                                                       *(hoelzlw), hoelzlw AT in.tum.de
                                                                                       *
                                                                                       **************************************************************************************/

#pragma once

#ifndef SIMD_TYPES_H
#define SIMD_TYPES_H

/*
 * Check, whether the function definitions are already included.
 * If yes, this denotes an error.
 *
 * The SIMD_DEFINITIONS.hpp-file needs the control macros set in this file to work properly.
 */
#ifdef SIMD_DEFINITIONS_H
#error \
    "SIMD Definitions already included! Never include that file directly! Include it only via including SIMD_TYPES (this file!)"
#endif /* defined SIMD_DEFINITIONS_H */

/*
 * Check, whether a solver is chosen, that uses the macros in this file.
 */
#if not WAVE_PROPAGATION_SOLVER == 5
#pragma message "SIMD macros included but non-vectorized solver specified"
#endif /* not WAVE_PROPAGATION_SOLVER == 5  */

/*
 * Care about precision.
 * Use single precision as default.
 *
 * Additionally, declare the macro SHIFT_SIGN_RIGHT to work properly with 32 and 64 bit.
 * Note the INTENTIONALLY FORGOTTEN semicolon at the end of the SHIFT_SIGN_RIGHT-definitions.
 * This forces the user to write the semicolon himself!
 */
#if defined FLOAT64
#pragma message "Using double as type for real numbers"
typedef double real;
#pragma message "Using unsigned long long as type for integer numbers"
typedef unsigned long long integer;

#define SHIFT_SIGN_RIGHT(x) \
    static_cast<integer>(static_cast<integer>(1) << (static_cast<integer>(64) - static_cast<integer>(x)))
#else /* not defined FLOAT64  */
#pragma message "Using float as type for real numbers"
typedef float        real;
#pragma message "Using unsigned int as type for integer numbers"
typedef unsigned int integer;

#define SHIFT_SIGN_RIGHT(x) (1 << (32 - static_cast<integer>(x)))
#endif /* not defined FLOAT64  */

/*
 * Set control macros
 *
 * Declare the vector length
 * Additionally declare a configuration specific macro of the form
 * 	VECTOR_extension_precision
 *
 * Moreover, declare an integer, representing the number,
 * which is returned by the instruction MOVEMASK, if the instruction is called on a vector with ALL SIGN BITS set
 * Use is as
 *
 * 	const real_vector vector = CONDITION(operand_a, operand_b);
 *
 * 	if (MOVEMASK(vector) == VECTOR_FULL_MASK) {
 *		// all components fullfill the condition
 *	} else {
 *		// some components do not fullfill the condition
 *	}
 */
#if (defined __SSE4_1__ and not defined __AVX__) or (defined __AVX__ and defined AVX128)
#pragma message "Using SSE4.1 for vectorization"
#include <smmintrin.h>

typedef __m128i integer_vector;
#if defined     FLOAT64
#define VECTOR_LENGTH 2
#define VECTOR_SSE4_FLOAT64
#define VECTOR_FULL_MASK 0x00000003

typedef __m128d real_vector;
#pragma message "Using vectors of 2 doubles"
#else /* not defined FLOAT64  */
#define VECTOR_LENGTH 4
#define VECTOR_SSE4_FLOAT32
#define VECTOR_FULL_MASK 0x0000000F

typedef __m128 real_vector;
#pragma message "Using vectors of 4 floats"
#endif /* not defined FLOAT64 */
#elif defined __AVX__
#pragma message "Using AVX for vectorization"
#include <immintrin.h>

typedef __m256i integer_vector;
#if defined FLOAT64
#define VECTOR_LENGTH 4
#define VECTOR_AVX_FLOAT64
#define VECTOR_FULL_MASK 0x0000000F

typedef __m256d real_vector;
#pragma message "Using vectors of 4 doubles"
#else /* not defined FLOAT64  */
#define VECTOR_LENGTH 8
#define VECTOR_AVX_FLOAT32
#define VECTOR_FULL_MASK 0x000000FF

typedef __m256 real_vector;
#pragma message "Using vectors of 8 floats"
#endif /* not defined FLOAT64  */
#else  /* not defined __SSE4__ and not defined __AVX__  */
#pragma message "Using no vectorization at all"
#define VECTOR_LENGTH 1
#define VECTOR_NOVEC
#endif /* not defined __SSE4__ and not defined __AVX__ */

/*
 * Control macros are set
 *
 * Include the function macros
 */
#include "SIMDDefinitions.hpp"

/*
 * Include the cost macros if flop counting is demanded
 */
#if defined COUNTFLOPS
#include "SIMD_COSTS.hpp"
#endif /* defined COUNTFLOPS  */

#endif /* #ifndef SIMD_TYPES_H */
