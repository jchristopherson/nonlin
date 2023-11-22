# If found, use the installed version; else, import the library
if (${LAPACK_FOUND})
    # Inform the user of what's going on
    message(STATUS "LAPACK library found.")
else()
    # Inform the user of what's going on
    message(STATUS "LAPACK not found.  Downloading the reference LAPACK.")

    # Fetch the proper content
    FetchContent_Declare(
        LAPACK
        GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
    )

    FetchContent_MakeAvailable(LAPACK)

    if (WIN32)
        if (${BUILD_SHARED_LIBS})
            add_custom_command(
                TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:lapack>
                $<TARGET_FILE_DIR:${PROJECT_NAME}
            )
        endif()
    endif()
endif()