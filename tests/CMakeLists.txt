include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
add_subdirectory(fortran_test_helper)

# Define the source files
set(nonlin_test_sources
    nonlin_test.f90
    nonlin_test_jacobian.f90
    nonlin_test_solve.f90
    nonlin_test_poly.f90
    nonlin_test_optimize.f90
    powell_badly_scaled.f90
)

add_executable(nonlin_test ${nonlin_test_sources})
link_library(nonlin_test ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(nonlin_test ${fortran_test_helper_LIBRARY} ${fortran_test_helper_INCLUDE_DIR})
add_test(
    NAME nonlin_test
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:nonlin_test>
)
