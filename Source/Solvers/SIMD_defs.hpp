#pragma once

#if defined(ENABLE_VECTORIZATION) && defined(ENABLE_VECTORIZATION_WITH_SIMD)
#define STRIDE 4u
#else
#define STRIDE 1u
#endif

#define PTR(var) (&var)