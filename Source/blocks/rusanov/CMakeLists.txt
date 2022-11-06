if( NOT ENABLE_CUDA )
if(SOLVER STREQUAL "rusanov")
    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT SWE_RusanovBlock.cpp)
    set(SELECTED_SWE_BLOCK_HEADER SWE_RusanovBlock.hh)
endif()
endif()

target_sources(swe 
    PRIVATE
    ${SELECTED_SWE_BLOCK_COMPILATION_UNIT}
    PUBLIC
    ${SELECTED_SWE_BLOCK_HEADER}
)


