#Todo: add other solvers

if (NOT ENABLE_CUDA)
    if (SOLVER STREQUAL "hlle")
        set(SELECTED_SWE_BLOCK_COMPILATION_UNIT SWE_WaveAccumulationBlock.cpp)
        set(SELECTED_SWE_BLOCK_HEADER SWE_WaveAccumulationBlock.hh)
    elseif (SOLVER STREQUAL "fwave")
        set(SELECTED_SWE_BLOCK_COMPILATION_UNIT SWE_WaveAccumulationBlock.cpp)
        set(SELECTED_SWE_BLOCK_HEADER SWE_WaveAccumulationBlock.hh)
    elseif(SOLVER STREQUAL "rusanov")
        set(SELECTED_SWE_BLOCK_COMPILATION_UNIT rusanov/SWE_RusanovBlock.cpp)
        set(SELECTED_SWE_BLOCK_HEADER rusanov/SWE_RusanovBlock.hh)
    elseif(SOLVER STREQUAL "augRie")
        set(SELECTED_SWE_BLOCK_COMPILATION_UNIT SWE_WaveAccumulationBlock.cpp)
        set(SELECTED_SWE_BLOCK_HEADER SWE_WaveAccumulationBlock.hh)
    else()
        message (WARNING "${SOLVER} not recognized!")
    endif()
endif()

target_sources(swe 
    PRIVATE
    SWE_Block.cpp
    ${SELECTED_SWE_BLOCK_COMPILATION_UNIT}
    PUBLIC
    SWE_Block.hh
    ${SELECTED_SWE_BLOCK_HEADER}
)
