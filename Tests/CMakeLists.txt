file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*.cpp")

foreach(file ${SOURCES})
  get_filename_component(filename ${file} NAME_WLE)
  add_executable(${filename} ${file})
  add_test(NAME ${filename} COMMAND ${filename})
  target_link_libraries(${filename} PRIVATE ${SWE_PROJECT_NAME} Catch2 Catch2WithMain)
endforeach()
