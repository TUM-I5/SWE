# include(FetchContent)

# FetchContent_Declare(
#     SWE-Solvers
#     GIT_REPOSITORY https://github.com/TUM-I5/SWE-Solvers.git
#     GIT_TAG        master
# )
# FetchContent_MakeAvailable(SWE-Solvers)

set(SOLVERS_LIB ${META_PROJECT_NAME}_SOLVERS)
add_library(${SOLVERS_LIB} INTERFACE)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "Solvers/*")

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})
target_sources(${META_PROJECT_NAME} INTERFACE ${SOURCES})
target_include_directories(${SOLVERS_LIB} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Solvers)

# target_include_directories(${SOLVERS_LIB} 
#     INTERFACE
#         "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers>"
#         "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

option(ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS "Enable augmented Riemann eigen coefficients" OFF)
if(ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS)
    target_compile_definitions(${SOLVERS_LIB} INTERFACE ENABLE_AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS)
endif()

# add_custom_target(${SOLVERS_LIB}- SOURCES ${SOURCES})
