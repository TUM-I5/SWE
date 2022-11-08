cmake_minimum_required(VERSION 3.20)

if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH_SHORT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
    set(GIT_BRANCH "")
    set(GIT_COMMIT_HASH "")
    set(GIT_COMMIT_HASH_SHORT "")
endif()

set(META_PROJECT_NAME        "SWE")
set(META_PROJECT_DESCRIPTION "An Education-Oriented Code for Parallel Tsunami Simulation")
set(META_AUTHOR_ORGANIZATION "Chair of Scientific Computing, Department of Informatics, Technische Universitaet Muenchen")
set(META_AUTHOR_DOMAIN       "http://www5.in.tum.de/SWE")
set(META_VERSION_REVISION    "${GIT_COMMIT_HASH_SHORT}")
set(META_GIT_BRANCH          "${GIT_BRANCH}")
set(META_GIT_HASH            "${GIT_COMMIT_HASH}")
message(STATUS "On Git Branch: ${GIT_BRANCH} (${GIT_COMMIT_HASH})")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_INLINES_HIDDEN YES)
set(CMAKE_CXX_EXTENSIONS OFF)
set(META_COMPILER_VERSION "${CMAKE_SYSTEM_NAME} ${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_COMPILER_VERSION})")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT ${META_PROJECT_NAME})

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

set(OPT_WARNINGS_AS_ERRORS OFF CACHE BOOL "" FORCE)
include(FetchContent)
FetchContent_Declare(
    _project_options
    GIT_REPOSITORY https://github.com/aminya/project_options.git
    GIT_TAG main
)
FetchContent_MakeAvailable(_project_options)

project(${META_PROJECT_NAME}
        DESCRIPTION ${META_PROJECT_DESCRIPTION}
        HOMEPAGE_URL ${META_AUTHOR_DOMAIN}
        LANGUAGES C CXX
)
file(WRITE "${PROJECT_BINARY_DIR}/.gitignore" "*")

include(${_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake)
dynamic_project_options()
target_compile_features(project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

include(GNUInstallDirs)
set(DEFAULT_CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

add_subdirectory(CMake)
include(MessageUtils)
include(Clang-Format)

include(TestBigEndian)
test_big_endian(BIG_ENDIAN)
include(ProcessorCount)
ProcessorCount(AVAILABLE_PROCESSOR_COUNT)
find_package(Threads REQUIRED)

enable_testing()
add_subdirectory(Source)
add_subdirectory(Tests)

display_separator("Project Configuration Summary")

display_variable(AVAILABLE_PROCESSOR_COUNT)

display_variable(CMAKE_BUILD_TYPE)
display_variable(CMAKE_CONFIGURATION_TYPES)
display_variable(CMAKE_CXX_FLAGS)
display_variable(CMAKE_CXX_FLAGS_DEBUG)
display_variable(CMAKE_CXX_FLAGS_RELEASE)
display_variable(CMAKE_CXX_FLAGS_MINSIZEREL)
display_variable(CMAKE_CXX_FLAGS_RELWITHDEBINFO)
display_variable(CMAKE_VERBOSE_MAKEFILE)

display_subseparator("Options to Build")
display_variable(BUILD_SHARED_LIBS)

display_subseparator("Options to Install")
display_variable(CMAKE_INSTALL_FULL_BINDIR)
display_variable(CMAKE_INSTALL_FULL_LIBDIR)

display_subseparator("End of Summary")
