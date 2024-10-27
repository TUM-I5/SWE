cmake_minimum_required(VERSION 3.21)

set(SWE_PROJECT_NAME        "SWE")
set(SWE_PROJECT_DESCRIPTION "An Education-Oriented Code for Tsunami Simulation")
set(SWE_AUTHOR_ORGANIZATION "Chair of Scientific Computing, Department of Informatics, Technische Universitaet Muenchen")
set(SWE_AUTHOR_DOMAIN       "http://www5.in.tum.de/SWE")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

include(GetGitRevisionDescription)
get_git_head_revision(GIT_BRANCH GIT_COMMIT_HASH)
get_git_commit_message(GIT_COMMIT_MESSAGE)
get_git_commit_timestamp(GIT_COMMIT_TIMESTAMP)

set(SWE_VERSION_REVISION "${GIT_COMMIT_HASH}")
set(SWE_GIT_BRANCH "${GIT_BRANCH}")
set(SWE_GIT_HASH "${GIT_COMMIT_HASH}")
set(SWE_GIT_COMMIT_MESSAGE "${GIT_COMMIT_MESSAGE}")
set(SWE_GIT_URL "https://gitlab.lrz.de/hpcsoftware/Peano")
message(STATUS "Revision: ${GIT_BRANCH} -m \"${GIT_COMMIT_MESSAGE}\" (${GIT_COMMIT_HASH}) ${GIT_COMMIT_TIMESTAMP}")

project(${SWE_PROJECT_NAME}
  DESCRIPTION ${SWE_PROJECT_DESCRIPTION}
  HOMEPAGE_URL ${SWE_AUTHOR_DOMAIN}
  LANGUAGES C CXX
)

set(CMAKE_C_STANDARD 17 CACHE STRING "")
set(CMAKE_C_STANDARD_REQUIRED TRUE)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "")
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_INLINES_HIDDEN YES)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON CACHE INTERNAL "")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

file(WRITE "${PROJECT_BINARY_DIR}/.gitignore" "*")

include(GNUInstallDirs)
set(DEFAULT_CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

include(TestBigEndian)
test_big_endian(BIG_ENDIAN)

include(ProcessorCount)
ProcessorCount(AVAILABLE_PROCESSOR_COUNT)

include(CTest)
enable_testing()

include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)
if(IPO_SUPPORTED)
  message(STATUS "IPO/LTO enabled")
  set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
else()
  #message(WARNING "IPO/LTO not supported: ${IPO_ERROR}")
  message(WARNING "IPO/LTO not supported")
endif()

add_library(SWE-Interface INTERFACE)
target_compile_features(SWE-Interface INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

target_compile_options(SWE-Interface INTERFACE
  -W
  -Wall
  -Wextra
  -Wpedantic
)

option(ENABLE_SINGLE_PRECISION "Enable single floating-point precision" OFF)
if(ENABLE_SINGLE_PRECISION)
  target_compile_definitions(SWE-Interface INTERFACE ENABLE_SINGLE_PRECISION)
endif()

option(ENABLE_OPENMP "Enable shared-memory parallelization of intra-block iterations using OpenMP" OFF)
if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(SWE-Interface INTERFACE OpenMP::OpenMP_CXX)
  target_compile_definitions(SWE-Interface INTERFACE ENABLE_OPENMP)
endif()

option(ENABLE_MPI "Enable distributed-memory parallelization using MPI" ON)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(SWE-Interface INTERFACE MPI::MPI_CXX)
  target_compile_definitions(SWE-Interface INTERFACE SWE_ENABLE_MPI)
  target_compile_definitions(SWE-Interface INTERFACE ENABLE_MPI)
endif()

option(ENABLE_GPI "Enable distributed-memory parallelization using GPI" OFF)
if(ENABLE_GPI)
  message(STATUS "GPI2 install path ${GPI2_INSTALL_PATH}")
  find_library(GPI2_LIB NAMES GPI2 REQUIRED HINTS ${GPI2_INSTALL_PATH}/lib64)
  find_path(GPI2_INC NAMES GASPI.h REQUIRED HINTS ${GPI2_INSTALL_PATH}/include)
  message(STATUS "GPI2 LIB INC: ${GPI2_LIB} ${GPI2_INC}")
  target_link_libraries(SWE-Interface INTERFACE ${GPI2_LIB})
  target_include_directories(SWE-Interface INTERFACE ${GPI2_INC})
  target_compile_definitions(SWE-Interface INTERFACE ENABLE_GPI)
endif()

option(ENABLE_VECTORIZATION "Enable auto-vectorization of the block iteration on the CPU (only supported solvers), and in Release Mode" OFF)
if(ENABLE_VECTORIZATION)
  target_compile_options(SWE-Interface INTERFACE
    $<$<CXX_COMPILER_ID:Clang>: -march=native -fno-math-errno>
    $<$<CXX_COMPILER_ID:AppleClang>: -march=native -fno-math-errno>
    $<$<CXX_COMPILER_ID:Intel>: -xHost>
    $<$<CXX_COMPILER_ID:GNU>: -march=native -ffast-math -fassociative-math>
  )
  target_compile_definitions(SWE-Interface INTERFACE ENABLE_VECTORIZATION)
endif()

option(ENABLE_NETCDF "Enable output writing with NetCDF." ON)
if(ENABLE_NETCDF)
  set(NETCDF_CXX "YES")
  find_package(NetCDF REQUIRED)
  target_link_libraries(SWE-Interface INTERFACE ${NETCDF_LIBRARIES} ${NETCDF_LIBRARIES_CXX})
  target_include_directories(SWE-Interface INTERFACE ${NETCDF_INCLUDES})
  target_compile_definitions(SWE-Interface INTERFACE ENABLE_NETCDF)
endif()

list(APPEND AVAILABLE_SOLVERS
  "FWave"
  "HLLE"
  "AugRie"
  "Rusanov"
)
set(WITH_SOLVER "FWave" CACHE STRING "Choose which Riemann solver to use (${AVAILABLE_SOLVERS})")
set_property(CACHE WITH_SOLVER PROPERTY STRINGS ${AVAILABLE_SOLVERS})

if(WITH_SOLVER STREQUAL "FWave")
  target_compile_definitions(SWE-Interface INTERFACE WITH_SOLVER_FWAVE)
elseif(WITH_SOLVER STREQUAL "HLLE")
  target_compile_definitions(SWE-Interface INTERFACE WITH_SOLVER_HLLE)
elseif(WITH_SOLVER STREQUAL "AugRie")
  target_compile_definitions(SWE-Interface INTERFACE WITH_SOLVER_AUGRIE)
elseif(WITH_SOLVER STREQUAL "Rusanov")
  target_compile_definitions(SWE-Interface INTERFACE WITH_SOLVER_RUSANOV)
endif()

find_package(Catch2 REQUIRED)
find_package(SWE-Solvers REQUIRED)

add_subdirectory(Source)
add_subdirectory(Tests)
