add_library(${SWE_PROJECT_NAME} OBJECT)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*")
list(FILTER SOURCES EXCLUDE REGEX ".*Runner\\.cpp$")

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})

target_sources(${SWE_PROJECT_NAME} PRIVATE ${SOURCES})

target_link_libraries(${SWE_PROJECT_NAME} PUBLIC SWE-Interface SWE-Solvers)
target_include_directories(${SWE_PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_MPI)
  add_executable(${SWE_PROJECT_NAME}-MPI-Runner Runners/MPI-Runner.cpp)
  target_link_libraries(${SWE_PROJECT_NAME}-MPI-Runner PRIVATE ${SWE_PROJECT_NAME})
endif()

if(ENABLE_GPI)
  add_executable(${SWE_PROJECT_NAME}-GPI-Runner Runners/GPI-Runner.cpp ${GPI2_INC}/GASPI.h)
  target_link_libraries(${SWE_PROJECT_NAME}-GPI-Runner PRIVATE ${SWE_PROJECT_NAME})
endif()
