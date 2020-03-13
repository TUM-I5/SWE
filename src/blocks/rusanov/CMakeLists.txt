if(SOLVER STREQUAL "rusanov")
    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT rusanov/SWE_RusanovBlock.cpp)
    set(SELECTED_SWE_BLOCK_HEADER rusanov/SWE_RusanovBlock.hh)
endif()

target_sources(swe 
    PRIVATE
    ${SELECTED_SWE_BLOCK_COMPILATION_UNIT}
    PUBLIC
    ${SELECTED_SWE_BLOCK_HEADER}
)


