#pragma once

#include <cuda_runtime.h>

#include "Tools/HelperFunctions.hpp"

/**
 * @brief Check the return value of the CUDA runtime API call and exit
 * the application if the call has failed.
 *
 * @param[in] msg Error message to print if the error code is not cudaSuccess
 * @param[in] ptr ptr for which the error is checked
 */
#define CUDA_CHECK_ERROR(msg, ptr) \
    { \
        cudaError_t err = cudaGetLastError(); \
        if (err == cudaError::cudaSuccess) \
        { \
            LOG_DBG(msg " " #ptr " ... passed!"); \
        } \
        else \
        { \
            LOG_FAT << msg " " #ptr " ... failed!\n[CUDA] " << cudaGetErrorString(err) << "\n"; \
        } \
    }

#ifdef ENABLE_CUDA_UMA_ALLOCATOR
// TODO - check if this is the correct way to initialize UMA memory
/**
 * @brief Allocate Unified Memory Access (UMA) memory on the device.
 *
 * @param[in] ptr Pointer to the memory to allocate
 * @param[in] size Size of the memory to allocate
 */
#define CUDA_MALLOC(ptr, size) \
    cudaMallocManaged(ptr, size); \
    CUDA_CHECK_ERROR("Allocating UMA memory ", #ptr);
#elif defined(ENABLE_CUDA_PINNEDMEM_ALLOCATOR)
#define CUDA_MALLOC(ptr, size) \
    cudaMallocHost(ptr, size); \
    CUDA_CHECK_ERROR("Allocating Pinned memory ", #ptr);
#else
/**
 * @brief Allocate memory on the device.
 *
 * @param[in] ptr Pointer to the memory to allocate
 * @param[in] size Size of the memory to allocate
 */
#define CUDA_MALLOC(ptr, size) \
    cudaMalloc(ptr, size); \
    CUDA_CHECK_ERROR("Allocating Device memory ", #ptr);
#endif

/**
 * @brief Copy from host/device memory to device/host.
 *
 * @param[in] to Pointer to the memory to copy to
 * @param[in] from Pointer to the memory to copy from
 * @param[in] size Size of the memory to copy
 * @param[in] direction Direction of the copy
 */
#define CUDA_MEM_CPY(to, from, size, direction, ...) \
    cudaMemcpy(to, from, size, direction); \
    CUDA_CHECK_ERROR("Copying " #from " to " #to " direction " #direction " ... " __VA_ARGS__, #to);

/**
 * @brief Freeing memory
 *
 * @param[in] ptr Pointer to the memory to free
 */
#define CUDA_FREE(ptr) \
    cudaFree(ptr); \
    CUDA_CHECK_ERROR("Freeing ", #ptr);

/**
 * @brief Check the return value of the CUDA runtime API call and exit
 * the application if the call has failed.
 *
 * @param[in] msg Error message to print if the error code is not cudaSuccess
 */
static void checkCUDAError(const char* msg)
{
    cudaError_t err = cudaGetLastError();
    if (cudaSuccess != err)
    {
        LOG_FAT << stderr, "\nCuda error (%s): %s.\n", msg, cudaGetErrorString(err);
        exit(-1);
    }
}

/**
 * @brief Check the return value of the CUDA runtime API call and exit
 * the application if the call has failed.
 *
 * @param[in] err Error code to check
 * @param[in] msg Error message to print if the error code is not cudaSuccess
 */
static void tryCUDA(cudaError_t err, const char* msg)
{
    if (cudaSuccess != err)
    {
        fprintf(stderr, "\nCuda error (%s): %s.\n", msg, cudaGetErrorString(err));
        exit(-1);
    }
}
