add_library(${META_PROJECT_NAME} OBJECT)
display_header("Creating Makefile of ${META_PROJECT_NAME}")

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*")
list(FILTER SOURCES EXCLUDE REGEX ".*Runner\\.cpp$")
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})
target_sources(${META_PROJECT_NAME} PRIVATE ${SOURCES})
target_include_directories(${META_PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

find_package(SWE-Solvers REQUIRED)
# target_link_libraries(${META_PROJECT_NAME}
#     PUBLIC
#         project_options
#         project_warnings
#         SWE-Solvers
# )

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
    include(CheckLanguage)
    check_language(CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_INCLUDE_DIRS "/usr/local/cuda/include")
    # uncomment for compiling outside singularity only
    # set(CMAKE_CUDA_INCLUDE_DIRS "/nfs/c3po/software/spack/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/cuda-11.4.4-5alqgf2yqocpxgtvqwqyzlywknrjte52/include")
    message(STATUS "CUDA enabled: ${CMAKE_CUDA_COMPILER_ID} ${CMAKE_CUDA_COMPILER} ${CMAKE_CUDA_COMPILER_VERSION} headers ${CMAKE_CUDA_INCLUDE_DIRS}")

    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    include_directories("${CMAKE_CUDA_INCLUDE_DIRS}")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_CUDA)
    # 86 = sm_86 = Ampere 2nd gen (RTX 30xx)
    # CUDA_STANDARD 17 is required for CUDA 11.0, CUDA20 fails to compile
    set_target_properties(${META_PROJECT_NAME} PROPERTIES CUDA_STANDARD 17 CUDA_ARCHITECTURES 86)
    option(CMAKE_CUDA_ARCHITECTURES 86) # for the cmake options

    # CUDA specific Options
    set(ENABLE_CUDA_ALLOCATOR "Device-Host" CACHE STRING "Choose which CUDA memory scheme to choose Device-Host, UMA, Pagelocked-PinnedMem")
    set(ENABLE_CUDA_ALLOCATOR_OPTIONS Device-Host UMA Pagelocked-PinnedMem)
    set_property(CACHE ENABLE_CUDA_ALLOCATOR PROPERTY STRINGS ${ENABLE_CUDA_ALLOCATOR_OPTIONS})

    if(ENABLE_CUDA_ALLOCATOR STREQUAL "UMA")
        message("CUDA UMA is chosen for memory.. Use CUDA Unified Memory Architecture Allocation instead of allocating data on Device.")
        target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_CUDA_UMA_ALLOCATOR)
    elseif(ENABLE_CUDA_ALLOCATOR STREQUAL "Pagelocked-PinnedMem")
        message("CUDA Page Locked memory is chosen for memory..")
        target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_CUDA_PINNEDMEM_ALLOCATOR)
    endif()

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
# if(OPT_ENABLE_SANITIZER_ADDRESS)
#     target_compile_options(${META_PROJECT_NAME} PRIVATE -fsanitize=address)
#     target_link_options(${META_PROJECT_NAME} PRIVATE -fsanitize=address)
# endif()

if(ENABLE_VECTORIZATION)
    target_compile_options(${META_PROJECT_NAME} PUBLIC
        $<$<CXX_COMPILER_ID:Clang>:  -march=core-avx2 -fno-math-errno>
        $<$<CXX_COMPILER_ID:AppleClang>: -march=native -fno-math-errno>
        $<$<CXX_COMPILER_ID:Intel>: -xHost -vec -simd> # TODO validate the flags to icpc
        $<$<CXX_COMPILER_ID:GNU>: -march=native -ffast-math -fassociative-math>
    )
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_VECTORIZATION)

    option(ENABLE_VECTORIZATION_WITH_SIMD "Use OMP vectorization or SIMD" ON)
    if(ENABLE_VECTORIZATION_WITH_SIMD)
        target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_VECTORIZATION_WITH_SIMD)
    endif()
endif()

option(ENABLE_WRITERS "Enable output writing" ON)
if(ENABLE_WRITERS)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_WRITERS)
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


option(ENABLE_TRACING_PROFILER "Enable Custom Profiling for SWE" ON)
if(ENABLE_TRACING_PROFILER)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_TRACING_PROFILER)
endif()

set(USED_SCENARIO "RadialDamBreakScenario" CACHE STRING "Choose a scenario to run")
set(USED_SCENARIO_OPTIONS BathymetryDamBreakScenario RadialDamBreakScenario SeaAtRestScenario SplashingConeScenario SplashingPoolScenario)
set_property(CACHE USED_SCENARIO PROPERTY STRINGS ${USED_SCENARIO_OPTIONS})
if(USED_SCENARIO STREQUAL "BathymetryDamBreakScenario")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC BATHYMETRYDAMBREAKSCENARIO)
elseif(USED_SCENARIO STREQUAL "RadialDamBreakScenario")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC RADIALDAMBREAKSCENARIO)
elseif(USED_SCENARIO STREQUAL "SeaAtRestScenario")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC SEAATRESTSCENARIO)
elseif(USED_SCENARIO STREQUAL "SplashingConeScenario")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC SPLASHINGCONESCENARIO)
elseif(USED_SCENARIO STREQUAL "SplashingPoolScenario")
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC SPLASHINGPOOLSCENARIO)
endif()

if(ENABLE_GPI)
    add_executable(${META_PROJECT_NAME}-GPI-Runner Runners/GPI-Runner.cpp ${GPI2_INC}/GASPI.h)
    target_link_libraries(${META_PROJECT_NAME}-GPI-Runner PRIVATE ${META_PROJECT_NAME})
    target_include_directories(${META_PROJECT_NAME}-GPI-Runner PRIVATE ${META_PROJECT_NAME})
endif()

