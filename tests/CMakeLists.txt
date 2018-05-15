# Define the source files
set(nonlin_test_sources
    ${PROJECT_SOURCE_DIR}/src/external/linalg/tests/test_core.f90
    nonlin_test.f90
    nonlin_test_jacobian.f90
    nonlin_test_solve.f90
    nonlin_test_poly.f90
    nonlin_test_optimize.f90
    powell_badly_scaled.f90
)

# Build the executable
add_executable(nonlin_test ${nonlin_test_sources})
target_link_libraries(nonlin_test nonlin)

# Run the tests
add_custom_command(
    OUTPUT nonlin_tests
    DEPENDS ${nonlin_test_sources}
    COMMAND nonlin_test
)
add_custom_target(run_nonlin_tests ALL DEPENDS nonlin_tests)

# Ensure LINALG, FERROR, and BLAS libraries are copied to the runtime directory
get_target_property(ferror_LibLocation ferror LOCATION)
add_custom_command(TARGET nonlin_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${ferror_LibLocation} $<TARGET_FILE_DIR:nonlin_test>
)

get_target_property(linalg_LibLocation linalg LOCATION)
add_custom_command(TARGET nonlin_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${linalg_LibLocation} $<TARGET_FILE_DIR:nonlin_test>
)

get_target_property(blas_LibLocation blas LOCATION)
add_custom_command(TARGET nonlin_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${blas_LibLocation} $<TARGET_FILE_DIR:nonlin_test>
)
