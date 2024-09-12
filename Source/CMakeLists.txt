add_library(${META_PROJECT_NAME} OBJECT)
display_header("Creating Makefile of ${META_PROJECT_NAME}")

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*")
list(FILTER SOURCES EXCLUDE REGEX ".*Runner\\.cpp$")
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})
target_sources(${META_PROJECT_NAME} PRIVATE ${SOURCES})
target_include_directories(${META_PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

find_package(SWE-Solvers REQUIRED)
target_link_libraries(${META_PROJECT_NAME}
    PUBLIC
        project_options
        project_warnings
        SWE-Solvers
)

option(ENABLE_SINGLE_PRECISION "Enable single floating-point precision" OFF)
if(ENABLE_SINGLE_PRECISION)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_SINGLE_PRECISION)
endif()

option(ENABLE_OPENMP "Enable shared-memory parallelization of intra-block iterations using OpenMP." OFF)
if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(${META_PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_OPENMP)
endif()

option(ENABLE_CUDA "Enable CUDA-accelerated SWE Blocks." OFF)
if(ENABLE_CUDA)
    message(FATAL_ERROR "CUDA Version not implemented")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_CUDA)
endif()

option(ENABLE_VISUALIZER "Enable the SWE-Visualizer which is built together with SWE." OFF)
if(ENABLE_VISUALIZER)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_VISUALIZER)
endif()

option(ENABLE_MPI "Enable distributed-memory parallelization using MPI." ON)
if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
    target_link_libraries(${META_PROJECT_NAME} PUBLIC MPI::MPI_CXX)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC SWE_ENABLE_MPI)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_MPI)
endif()

option(ENABLE_GPI "Enable distributed-memory parallelization using GPI." OFF)
if(ENABLE_GPI)
    message(STATUS "GPI2 install path ${GPI2_INSTALL_PATH}")
    find_library(GPI2_LIB NAMES GPI2 REQUIRED HINTS ${GPI2_INSTALL_PATH}/lib64)
    find_path(GPI2_INC NAMES GASPI.h REQUIRED HINTS ${GPI2_INSTALL_PATH}/include)
    message(STATUS "GPI2 LIB INC: ${GPI2_LIB} ${GPI2_INC}")
    target_link_libraries(${META_PROJECT_NAME} PUBLIC ${GPI2_LIB})
    target_include_directories(${META_PROJECT_NAME} PUBLIC ${GPI2_INC})
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC SWE_ENABLE_GPI)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_GPI)
endif()

option(ENABLE_VECTORIZATION "Enable auto-vectorization of the block iteration on the CPU (only supported solvers), and in Release Mode." OFF)
if(ENABLE_VECTORIZATION)
    target_compile_options(${META_PROJECT_NAME} PUBLIC
        $<$<CXX_COMPILER_ID:Clang>: -march=native -fno-math-errno>
        $<$<CXX_COMPILER_ID:AppleClang>: -march=native -fno-math-errno>
        $<$<CXX_COMPILER_ID:Intel>: -xHost>
        $<$<CXX_COMPILER_ID:GNU>: -march=native -ffast-math -fassociative-math>
    )
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_VECTORIZATION)
endif()

option(ENABLE_NETCDF "Enable output writing with NetCDF." ON)
if(ENABLE_NETCDF)
    if(NOT MSVC)
        set(NETCDF_CXX "YES")
    endif()
    find_package(NetCDF REQUIRED)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_NETCDF)
    #target_link_system_libraries(${META_PROJECT_NAME} PUBLIC ${NETCDF_LIBRARIES})
    if(NOT MSVC)
        target_link_libraries(${META_PROJECT_NAME} PUBLIC ${NETCDF_LIBRARIES} ${NETCDF_LIBRARIES_CXX})
    else()
        target_link_libraries(${META_PROJECT_NAME} PUBLIC ${NETCDF_LIBRARIES})
    endif()
    #target_include_system_directories(${META_PROJECT_NAME} PUBLIC ${NETCDF_INCLUDES})
    target_include_directories(${META_PROJECT_NAME} PUBLIC ${NETCDF_INCLUDES} )
endif()

set(WITH_SOLVER "FWave" CACHE STRING "Choose which Riemann solver to use (Default=FWave, AugRie, HLLE, Rusanov)")
set(WITH_SOLVER_OPTIONS FWave HLLE AugRie Rusanov)
set_property(CACHE WITH_SOLVER PROPERTY STRINGS ${WITH_SOLVER_OPTIONS})

if(WITH_SOLVER STREQUAL "FWave")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC WITH_SOLVER_FWAVE)
elseif(WITH_SOLVER STREQUAL "HLLE")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC WITH_SOLVER_HLLE)
elseif(WITH_SOLVER STREQUAL "AugRie")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC WITH_SOLVER_AUGRIE)
elseif(WITH_SOLVER STREQUAL "Rusanov")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC WITH_SOLVER_RUSANOV)
endif()

if(ENABLE_MPI)
    add_executable(${META_PROJECT_NAME}-MPI-Runner Runners/MPI-Runner.cpp)
    target_link_libraries(${META_PROJECT_NAME}-MPI-Runner PRIVATE ${META_PROJECT_NAME})
endif()

if(ENABLE_GPI)
    add_executable(${META_PROJECT_NAME}-GPI-Runner Runners/GPI-Runner.cpp ${GPI2_INC}/GASPI.h)
    target_link_libraries(${META_PROJECT_NAME}-GPI-Runner PRIVATE ${META_PROJECT_NAME})
    target_include_directories(${META_PROJECT_NAME}-GPI-Runner PRIVATE ${META_PROJECT_NAME})
endif()
