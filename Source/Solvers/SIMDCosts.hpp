/***********************************************************************************/ /**
                                                                                       *
                                                                                       * \file SIMD_COSTS.hpp
                                                                                       *
                                                                                       * \brief Contains latency costs
                                                                                       *for the intrinsics used for the
                                                                                       *vectorization
                                                                                       *
                                                                                       * \author Wolfgang Hölzl
                                                                                       *(hoelzlw), hoelzlw AT in.tum.de
                                                                                       *
                                                                                       **************************************************************************************/

#pragma once

#ifndef SIMD_COSTS_H
#define SIMD_COSTS_H

#ifdef COUNTFLOPS

#ifndef SIMD_TYPES_H
#error "SIMD_COSTS included without SIMD_TYPES! Never include this file directly! Include it only via SIMD_TYPES!"
#endif /* defined SIMD_TYPES_H */

#define COSTS_ADDS 1
#define COSTS_SUBS 1
#define COSTS_MULS 1
#define COSTS_DIVS 1
#define COSTS_SQRTS 1

#define COSTS_MAXS 1
#define COSTS_MINS 1

#define COSTS_CMPS 0

#define COSTS_ANDS 0
#define COSTS_ORS 0

#define COSTS_FABSS 1

#if defined VECTOR_SSE4_FLOAT32
#define COSTS_ADDV 4
#define COSTS_SUBV 4
#define COSTS_MULV 4
#define COSTS_DIVV 4
#define COSTS_SQRTV 4

#define COSTS_LOADU 0
#define COSTS_STOREU 0
#define COSTS_SETV_R 0
#define COSTS_SETV_I 0
#define COSTS_ZEROV_R 0
#define COSTS_ZEROV_I 0

#define COSTS_MAXV 4
#define COSTS_MINV 4

#define COSTS_CMP_LT 0
#define COSTS_CMP_LE 0
#define COSTS_CMP_GT 0
#define COSTS_CMP_GE 0

#define COSTS_CMP_EQ_I 0

#define COSTS_ANDV_R 0
#define COSTS_ORV_R 0
#define COSTS_XORV_R 0
#define COSTS_ORV_I 0
#define COSTS_ANDNOTV_R 0

#define COSTS_BLENDV 0
#define COSTS_BLENDV_I 0
#define COSTS_MOVEMASK 0
#define COSTS_SHIFT_LEFT 0

#define COSTS_FABS 4
#define COSTS_NOTV_R 0
#define COSTS_NOTV_I 0
#elif defined VECTOR_SSE4_FLOAT64
#error "SSE4 with double precision not implemented at the moment"
#elif defined VECTOR_AVX_FLOAT32
#define COSTS_ADDV 8
#define COSTS_SUBV 8
#define COSTS_MULV 8
#define COSTS_DIVV 8
#define COSTS_SQRTV 8

#define COSTS_LOADU 0
#define COSTS_STOREU 0
#define COSTS_SETV_R 0
#define COSTS_SETV_I 0
#define COSTS_ZEROV_R 0
#define COSTS_ZEROV_I 0

#define COSTS_MAXV 8
#define COSTS_MINV 8

#define COSTS_CMP_LT 0
#define COSTS_CMP_LE 0
#define COSTS_CMP_GT 0
#define COSTS_CMP_GE 0
#define COSTS_CMP_EQ_I 0

#define COSTS_ANDV_R 0
#define COSTS_ORV_R 0
#define COSTS_XORV_R 0
#define COSTS_ORV_I 0
#define COSTS_ANDNOTV_R 0

#define COSTS_BLENDV 0
#define COSTS_BLENDV_I 0
#define COSTS_MOVEMASK 0
#define COSTS_SHIFT_LEFT 0

#define COSTS_FABS 8
#define COSTS_NOTV_R 0
#define COSTS_NOTV_I 0
#elif defined VECTOR_AVX_FLOAT64
#error "AVX with double precision not implemented at the moment"
#else /* no vectorization type defined */
#pragma message "SIMD-Costs included, but no Vector-Type defined."
#endif

#endif /* not defined COUNTFLOPS */

#endif /* #ifndef SIMD_COSTS_H */
