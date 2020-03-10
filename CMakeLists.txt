cmake_minimum_required(VERSION 3.13)

project(SWE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(ENABLE_NETCDF "Enable Output using NetCDF" ON)
option(ENABLE_OPENMP "Enable shared-memory parallelization of intra-block iterations using OpenMP" ON)

find_package(NetCDF)
find_package(OpenMP)
find_package(MPI)

add_executable(swe-single-node)
add_library(swe)

add_subdirectory(src/blocks)
add_subdirectory(src/scenarios)
add_subdirectory(src/tools)
add_subdirectory(src/writer)

target_sources(swe-single-node 
    PRIVATE
    src/examples/swe_simple.cpp
)

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
endif()




target_include_directories(swe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(swe-single-node PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_definitions(swe PUBLIC WAVE_PROPAGATION_SOLVER=4)

target_link_libraries(swe-single-node PUBLIC swe)

# Warnings
target_compile_options(swe PUBLIC -Wall -Wextra -pedantic)

#if(NOT BUILD_USING_UPCXX_WRAPPER)
#    target_include_directories(actorlib PUBLIC ${UPCXX_INCLUDES})
#    target_compile_definitions(actorlib PUBLIC ${UPCXX_DEFINES})
#endif()
#
#if (ACTORLIB_USE_OPENMP_TASKS) 
#    target_compile_definitions(actorlib PRIVATE ACTORLIB_USE_OPENMP_TASKS)
#    target_link_libraries(actorlib PRIVATE OpenMP::OpenMP_CXX)
#    message(STATUS "Actorlib Parallelization set to OpenMP Tasks")
#endif()
#
##target_compile_definitions(actorlib PRIVATE ACTORLIB_GIT_REVISION=${GIT_SHA1})
##target_compile_definitions(actorlib PRIVATE ACTORLIB_GIT_DATE=${GIT_DATE})
##target_compile_definitions(actorlib PRIVATE ACTORLIB_UPCXX_GASNET_CONDUIT=${UPCXX_GASNET_CONDUIT_NAME})
##target_compile_definitions(actorlib PRIVATE ACTORLIB_UPCXX_INSTALLATION=${UPCXX_GASNET_INSTALLATION_PATH})
##target_compile_definitions(actorlib PRIVATE ACTORLIB_GIT_COMMIT_MSG=${GIT_COMMIT_SUBJECT})
#
#if (ENABLE_MEMORY_SANITATION AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
#    message(STATUS "Memory sanitation enabled. Linking appropriate libraries.")
#    target_compile_options(actorlib PUBLIC -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -fsanitize-address-use-after-scope -Wuninitialized -Wall -Wextra -pedantic -g)
#    target_link_libraries(actorlib PUBLIC -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -g)
#    # For some old Ubuntu machines, you may need to use a different linker. Set the following variable: CMAKE_EXE_LINKER_FLAGS=-fuse-ld=gold
#elseif(ENABLE_MEMORY_SANITATION)
#    message(STATUS "Full Memory sanitation only supported on GNU Environment. Trying to use intel fallback.")
#    target_compile_options(actorlib PUBLIC -g -traceback -check-uninit)
#    target_link_libraries(actorlib PUBLIC -g -traceback -check-uninit)
#    
##endif()

##target_include_directories(task_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
#target_compile_definitions(task_test PRIVATE USE_OPENMP)
#
#target_include_directories(channel_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
#target_include_directories(square_root_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
##target_include_directories(pingpong_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

#target_include_directories(pond PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
##target_include_directories(pond PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/applications/pond)
#target_compile_definitions(pond PRIVATE WAVE_PROPAGATION_SOLVER=3)
#target_compile_options(pond PRIVATE -fstrict-aliasing -fargument-noalias)

##if(BUILD_RELEASE)
#message(STATUS "Configure Pond for Release mode")
    #    if(NOT IS_CROSS_COMPILING)
    #    target_compile_options(pond PRIVATE -march=native)
    #else()
    #        target_compile_options(pond PRIVATE -xMIC-AVX512)
    #    endif()
    #target_compile_options(pond PRIVATE -O3)
    #target_compile_definitions(pond PRIVATE NDEBUG VECTORIZE)
    #if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    #    target_compile_options(pond PRIVATE -qopt-report-phase=vec -qopt-report=2)
    #endif()
    #else()
    #message(STATUS "Configure Pond for Debug mode")
    #endif()

    #if (NOT BUILD_USING_UPCXX_WRAPPER) 
    #  foreach(LIBNAME ${UPCXX_LIBS})
    #    target_link_libraries(actorlib PUBLIC ${UPCXX_LIB_${LIBNAME}})
    #endforeach()
    #target_link_libraries(actorlib PUBLIC pthread)
    #else()
    #add_custom_command(TARGET actorlib PRE_BUILD COMMAND export UPCXX_CODEMODE=debug COMMAND export UPCXX_THREADMODE=par)
    #target_link_libraries(actorlib PUBLIC pthread)
    #endif()

    #if(NETCDF_FOUND)
    #message(STATUS "Enabled NetCDF IO in Pond")
    #target_compile_definitions(pond PRIVATE WRITENETCDF)
    #target_include_directories(pond PRIVATE ${NETCDF_INCLUDE_DIRS})
    #target_link_libraries (pond PRIVATE ${NETCDF_LIBRARIES})
    #else()
    #message(WARNING "Enabled VTK IO (Slow, fall-back, as NetCDF is unavailable)")
    #target_compile_definitions(pond PRIVATE WRITEVTK)
    #endif()

    #if (METIS_FOUND AND ENABLE_METIS_GRAPH_PARTITIONING)
    #target_include_directories(pond PRIVATE ${METIS_INCLUDE_DIRS})
    #target_link_libraries(pond PRIVATE ${METIS_LIBRARIES})
    #target_compile_definitions(pond PRIVATE METIS_PARTITIONING)
    #message(STATUS "Enable advanced Graph Partitioning using the Metis graph partitioning library.")
    #elseif (EIGEN3_FOUND AND ENABLE_ADVANCED_GRAPH_PARTITIONING)
    #target_include_directories(pond PRIVATE ${EIGEN3_INCLUDE_DIR})
    #target_compile_definitions(pond PRIVATE EIGEN3_PARTITIONING)
    #message(STATUS "Enable advanced Graph Partitioning using Eigen3 library.")
    #elseif(ENABLE_ADVANCED_GRAPH_PARTITIONING)
    #message(ERROR "Unable to enable advanced graph partitioning, as Eigen3 library cannot be found.")
    #else()
    #    message(STATUS "Use fallback method of simple graph partitioning.")
    #endif()

#if(NOT ENABLE_FILE_OUTPUT)
    #message(STATUS "Disable writing of time step files in Pond")
    #target_compile_definitions(pond PRIVATE NOWRITE)
    #endif()
#
#target_link_libraries(task_test PRIVATE OpenMP::OpenMP_CXX)
#target_link_libraries(channel_test PUBLIC actorlib)
#target_link_libraries(square_root_example PUBLIC actorlib)
#target_link_libraries(pingpong_example PUBLIC actorlib)

#target_link_libraries(pond PRIVATE OpenMP::OpenMP_CXX)
#target_link_libraries(pond PRIVATE actorlib)
#
#set_property(TARGET actorlib PROPERTY POSITION_INDEPENDENT_CODE ON)
#set_property(TARGET actorlib PROPERTY CXX_STANDARD 14)

#set_property(TARGET pond PROPERTY CXX_STANDARD 14)
#set_property(TARGET channel_test PROPERTY CXX_STANDARD 14)
#set_property(TARGET task_test PROPERTY CXX_STANDARD 14)
#set_property(TARGET pingpong_example PROPERTY CXX_STANDARD 14)
#set_property(TARGET square_root_example PROPERTY CXX_STANDARD 14)

# Export compiler commands for use in YouCompleteMe
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
