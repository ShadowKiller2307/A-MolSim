option(BUILD_DOXY "Create Doxygen Documentation" ON)

if (BUILD_DOXY)
    find_package(Doxygen)

    if (DOXYGEN_FOUND)
        add_custom_target(
                doc_doxygen
                COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
				
        )

        set_target_properties(doc_doxygen PROPERTIES EXCLUDE_FROM_ALL TRUE)
    endif (DOXYGEN_FOUND)
endif (BUILD_DOXY)