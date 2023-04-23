set(DOXYGEN_CALLER_GRAPH YES)
set(DOXYGEN_CALL_GRAPH YES)
set(DOXYGEN_EXTRACT_ALL YES)
find_package(Doxygen
    #REQUIRED dot
    OPTIONAL_COMPONENTS dot mscgen dia
)

if(Doxygen_FOUND)
    if (EXISTS ${CMAKE_SOURCE_DIR}/Documentation/Doxyfile)
        set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/Documentation/Doxyfile)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Documentation/Doxyfile)

        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

        if(NOT TARGET doxygen-docs)
            add_custom_target(doxygen-docs
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Documentation
                COMMENT "Generating Doxygen documentation."
                VERBATIM
            )
        endif()
    endif()
endif()
