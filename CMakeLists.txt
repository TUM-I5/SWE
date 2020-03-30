cmake_minimum_required(VERSION 3.13)

project(SWE LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(ENABLE_CUDA                  "Enable CUDA-accelerated SWE_Blocks" OFF)
option(ENABLE_NETCDF                "Enable Output using NetCDF" ON)
option(ENABLE_OPENMP                "Enable shared-memory parallelization of intra-block iterations using OpenMP" ON)
option(ENABLE_VECTORIZATION         "Enable auto-vectorization of the block iteration on the CPU (only supported solvers), and in Release Mode." OFF)
option(ENABLE_VECTORIZATION_VERBOSE "Enable verbose auto-vectorization reporting." OFF)

#todo: adapt these to supported solvers
set(SOLVER_STRINGS "rusanov" "augRie" "fwave" "hlle")

set(SOLVER "augRie" CACHE STRING "Default Riemann solver")
set_property(CACHE SOLVER PROPERTY STRINGS ${SOLVER_STRINGS})

if (ENABLE_CUDA)
    enable_language(CUDA)
endif()

if( NOT SOLVER IN_LIST SOLVER_STRINGS)
    message(FATAL_ERROR "solver must be one of ${SOLVER_STRINGS}")
endif()

message(STATUS "Using riemann solver ${SOLVER}")

find_package(NetCDF)
find_package(OpenMP)
find_package(MPI)

add_executable(swe-single-node)
add_library(swe)

if(SOLVER STREQUAL "augRie")
    target_compile_definitions(swe PUBLIC SOLVER_AUGRIE)
elseif(SOLVER STREQUAL "hlle")
    target_compile_definitions(swe PUBLIC SOLVER_HLLE)
elseif(SOLVER STREQUAL "fwave")
    target_compile_definitions(swe PUBLIC SOLVER_FWAVE)
elseif(SOLVER STREQUAL "rusanov")
    target_compile_definitions(swe PUBLIC SOLVER_RUSANOV)
endif()

add_subdirectory(src/blocks)
add_subdirectory(src/blocks/rusanov)
add_subdirectory(src/scenarios)
add_subdirectory(src/tools)
add_subdirectory(src/writer)
add_subdirectory(src/solvers)

if (ENABLE_CUDA)
    target_compile_definitions(swe PUBLIC CUDA)
    add_subdirectory(src/blocks/cuda)
endif()

target_sources(swe-single-node 
    PRIVATE
    src/examples/swe_simple.cpp
)
    
target_link_libraries(swe-single-node PUBLIC swe)

if (MPI_FOUND)
    add_executable(swe-mpi)
    target_sources(swe-mpi 
        PRIVATE
        src/examples/swe_mpi.cpp
    )
    target_link_libraries(swe-mpi PUBLIC swe)
    target_link_libraries(swe-mpi PRIVATE ${MPI_CXX_LIBRARIES})
    target_include_directories(swe-mpi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_include_directories(swe-mpi PRIVATE ${MPI_CXX_INCLUDE_DIRS})
endif()

if (OpenMP_FOUND AND ENABLE_OPENMP)
    target_compile_definitions(swe PUBLIC LOOP_OPENMP)
    target_link_libraries(swe PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP supported loop parallelization is enabled.")
else()
    message(STATUS "OpenMP loop parallelization is disabled.")
endif()

if (OpenMP_FOUND AND ENABLE_VECTORIZATION)
    #target_compile_definitions(swe PUBLIC $<$<CONFIG:Release>:VECTORIZE NDEBUG>)
    target_compile_definitions(swe PUBLIC VECTORIZE NDEBUG)
    target_link_libraries(swe PRIVATE OpenMP::OpenMP_CXX)
    target_compile_options(swe PRIVATE 
        $<$<CXX_COMPILER_ID:Clang>: -march=native -fno-math-errno>
        $<$<CXX_COMPILER_ID:AppleClang>: -march=native -fno-math-errno>
        $<$<CXX_COMPILER_ID:Intel>: -xHost>
        $<$<CXX_COMPILER_ID:GNU>: -march=native -ffast-math -fassociative-math>
    )
    if(ENABLE_VECTORIZATION_VERBOSE)
        target_compile_options(swe PRIVATE
            $<$<CXX_COMPILER_ID:AppleClang>: -Rpass-analysis=loop-vectorize>
            $<$<CXX_COMPILER_ID:Clang>: -Rpass-analysis=loop-vectorize>
            $<$<CXX_COMPILER_ID:Intel>: -qopt-report-phase=vec -qopt-report=2>
            $<$<CXX_COMPILER_ID:GNU>: -ftree-vectorizer-verbose=4 -fopt-info-all-vec>
        )
    endif()
    message(STATUS "OpenMP vectorization is enabled (Only for Release Mode).")

elseif(OpenMP_NOTFOUND AND ENABLE_VECTORIZATION)
    message(FATAL_ERROR "OpenMP is required to use Vectorization.")
else()
    message(STATUS "OpenMP vectorization is disabled.")
endif()

target_include_directories(swe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(swe-single-node PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

if (ENABLE_NETCDF AND NETCDF_FOUND)
    message(STATUS "Enabled NetCDF writer")
    target_compile_definitions(swe PRIVATE WRITENETCDF)
    target_include_directories(swe PRIVATE ${NETCDF_INCLUDE_DIRS})
    target_link_libraries(swe PRIVATE ${NETCDF_LIBRARIES})
else ()
    message(WARNING "Enabled VTK IO (Slow, fall-back, as NetCDF is unavailable)")
endif()


# Warnings

target_compile_options(swe PUBLIC 
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    $<$<COMPILE_LANGUAGE:CUDA>:-g>
)

#target_compile_options(swe PUBLIC <$<$ >:-g>) # -Wextra -pedantic)

#C++ Standard
set_property(TARGET swe PROPERTY CXX_STANDARD 11)
set_property(TARGET swe PROPERTY CUDA_STANDARD 11)
set_property(TARGET swe-mpi PROPERTY CXX_STANDARD 11)
set_property(TARGET swe-single-node PROPERTY CXX_STANDARD 11)

# Export compiler commands for use in YouCompleteMe
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
