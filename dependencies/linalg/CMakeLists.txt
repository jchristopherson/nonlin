# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    linalg
    GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
)

FetchContent_MakeAvailable(linalg)

if (WIN32)
    if (BUILD_SHARED_LIBS)
        add_custom_command(
            TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:linalg>
            $<TARGET_FILE_DIR:${PROJECT_NAME}
        )
    endif()
endif()

set(linalg_INCLUDE_DIR ${linalg_BINARY_DIR}/include)
set(linalg_INCLUDE_DIR ${linalg_INCLUDE_DIR} PARENT_SCOPE)
configure_file(
    "${linalg_SOURCE_DIR}/include/linalg.h"
    "${linalg_INCLUDE_DIR}/linalg.h"
    COPYONLY
)

# Make a parent-scope variable for the library
set(linalg_LIBRARY linalg)
set(linalg_LIBRARY ${linalg_LIBRARY} PARENT_SCOPE)
