# Define the source files
set(nonlin_test_sources
    nonlin_test.f90
    nonlin_test_jacobian.f90
    nonlin_test_solve.f90
    nonlin_test_poly.f90
    nonlin_test_optimize.f90
    nonlin_test_simple.f90
    nonlin_test_iteration.f90
    nonlin_test_continuation.f90
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
    set(_linalg_build_testing ${BUILD_TESTING})
    set(BUILD_TESTING OFF)
    FetchContent_MakeAvailable(fortran_test_helper)
    set(BUILD_TESTING ${_linalg_build_testing})
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
set_tests_properties(nonlin_test PROPERTIES TIMEOUT 30)
add_test(
    NAME nonlin_sparse_jacobian
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:nonlin_test> sparse_jacobian
)
set_tests_properties(nonlin_sparse_jacobian PROPERTIES TIMEOUT 30)
foreach(mode enabled disabled)
    add_test(
        NAME nonlin_sparse_missing_${mode}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND ${CMAKE_COMMAND}
            "-DTEST_EXECUTABLE=$<TARGET_FILE:nonlin_test>"
            "-DTEST_ARGUMENT=sparse_missing_${mode}"
            "-DEXPECTED_ERROR=211"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/check_error_stop.cmake"
    )
    set_tests_properties(nonlin_sparse_missing_${mode} PROPERTIES TIMEOUT 30)
endforeach()
add_test(
    NAME nonlin_newton_line_search
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:nonlin_test> newton_line_search
)
set_tests_properties(nonlin_newton_line_search PROPERTIES TIMEOUT 30)
foreach(mode invalid_scale nan_scale invalid_fold zero_direction)
    add_test(
        NAME nonlin_newton_${mode}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND ${CMAKE_COMMAND}
            "-DTEST_EXECUTABLE=$<TARGET_FILE:nonlin_test>"
            "-DTEST_ARGUMENT=newton_${mode}"
            "-DEXPECTED_ERROR=201"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/check_error_stop.cmake"
    )
    set_tests_properties(nonlin_newton_${mode} PROPERTIES TIMEOUT 30)
endforeach()
add_test(
    NAME nonlin_jacobian_threads
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:nonlin_test> jacobian_threads
)
set_tests_properties(nonlin_jacobian_threads PROPERTIES
    TIMEOUT 30
    ENVIRONMENT "OMP_DYNAMIC=FALSE;OMP_THREAD_LIMIT=8"
)

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET nonlin_test
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:nonlin_test> 
        $<TARGET_FILE_DIR:nonlin_test>
        COMMAND_EXPAND_LISTS
    )
endif()
