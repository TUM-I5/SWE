#Todo: add other solvers


if (SOLVER STREQUAL "fwave")
    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA SWE_WavePropagationBlockCuda.cu)
    set(SELECTED_SWE_BLOCK_HEADER_CUDA SWE_WavePropagationBlockCuda.hh)

    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA_KERNELS SWE_WavePropagationBlockCuda_kernels.cu)
    set(SELECTED_SWE_BLOCK_HEADER_CUDA_KERNELS SWE_WavePropagationBlockCuda_kernels.hh)
elseif(SOLVER STREQUAL "rusanov")
    #Todo
elseif(SOLVER STREQUAL "augRie")
    #Todo
else()
    message (FATAL_ERROR "${SOLVER} is not supported with cuda")
endif()

target_sources(swe 
    PRIVATE
    SWE_BlockCUDA.cu
    SWE_BlockCUDA_kernels.cu
    ${SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA}
    ${SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA_KERNELS}
    PUBLIC
    SWE_BlockCUDA.hh
    SWE_BlockCUDA_kernels.hh
    ${SELECTED_SWE_BLOCK_HEADER_CUDA}
    ${SELECTED_SWE_BLOCK_HEADER_CUDA_KERNELS}
)
