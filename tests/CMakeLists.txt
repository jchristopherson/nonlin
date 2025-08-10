# Define the source files
set(nonlin_test_sources
    nonlin_test.f90
    nonlin_test_jacobian.f90
    nonlin_test_solve.f90
    nonlin_test_poly.f90
    nonlin_test_optimize.f90
    powell_badly_scaled.f90
)

# Dependencies
find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    FetchContent_MakeAvailable(fortran_test_helper)
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

add_executable(nonlin_test ${nonlin_test_sources})
target_link_libraries(
    nonlin_test
    nonlin
    ${FTH_LIBRARIES}
)
add_test(
    NAME nonlin_test
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:nonlin_test>
)

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET linalg_test
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:linalg_test> 
        $<TARGET_FILE_DIR:linalg_test>
        COMMAND_EXPAND_LISTS
    )
endif()
