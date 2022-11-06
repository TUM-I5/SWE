find_program(CLANGFORMAT clang-format NAMES clang-format)
if(CLANGFORMAT)
    execute_process(COMMAND ${CLANGFORMAT} --version
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        RESULT_VARIABLE CLANGFORMATRESULT
        OUTPUT_VARIABLE CLANGFORMATVERSION
        ERROR_VARIABLE CLANGFORMATERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(CLANGFORMATRESULT EQUAL 0)
        message(STATUS "Found ${CLANGFORMATVERSION}: " ${CLANGFORMAT})

        file(GLOB_RECURSE ALL_CXX_SOURCE_FILES
            ${CMAKE_SOURCE_DIR}/Source/*.[ch]pp
            ${CMAKE_SOURCE_DIR}/Source/*.cpph
            ${CMAKE_SOURCE_DIR}/Tests/*.[ch]pp
            ${CMAKE_SOURCE_DIR}/Tests/*.cpph
        )

        if(NOT TARGET clang-format)
            add_custom_target(clang-format
                COMMAND ${CLANGFORMAT}
                    -i
                    -style=file
                    -fallback-style=none
                    -verbose
                    ${ALL_CXX_SOURCE_FILES}
                SOURCES "${CMAKE_SOURCE_DIR}/.clang-format"
                COMMENT "Format all source files. This may take a while..."
            )
        endif()

        if(NOT TARGET clang-format-check)
            add_custom_target(clang-format-check
                # Use ! to negate the result for correct output
                COMMAND !
                    ${CLANGFORMAT}
                    -style=file
                    -output-replacements-xml
                    -fallback-style=none
                    -verbose
                    ${ALL_CXX_SOURCE_FILES}
                    | grep -q "Replacement offset"
                SOURCES "${CMAKE_SOURCE_DIR}/.clang-format"
                COMMENT "Checking clang-format changes."
            )
        endif()

        if(NOT TARGET clang-format-dry)
            add_custom_target(clang-format-dry
                COMMAND ${CLANGFORMAT}
                    -style=file
                    -dry-run
                    -fallback-style=none
                    ${ALL_CXX_SOURCE_FILES}
                SOURCES "${CMAKE_SOURCE_DIR}/.clang-format"
                COMMENT "Running clang-format in dry mode."
            )
        endif()
    else()
        message(SEND_ERROR "clang-format found but cannot retrieve version information")
    endif()
else()
    message(STATUS "clang-format requested but executable not found")
endif()
