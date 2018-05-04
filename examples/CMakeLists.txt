# Quasi-Newton Example
add_executable(quasi_newton_example nonlin_quasi_newton_example.f90)
target_link_libraries(quasi_newton_example nonlin)

# Least-Squares Example
add_executable(least_squares_example nonlin_least_squares_example.f90)
target_link_libraries(least_squares_example nonlin)

# Polynomial Fitting Example
add_executable(poly_fit_example nonlin_poly_fit_example.f90)
target_link_libraries(poly_fit_example nonlin)

# Nelder-Mead Example
add_executable(nm_opt_example nonlin_nm_opt_example.f90)
target_link_libraries(nm_opt_example nonlin)

# Newton Example
add_executable(newton_example nonlin_newton_example.f90)
target_link_libraries(newton_example nonlin)


# Ensure LINALG, FERROR, and BLAS libraries are copied to the runtime directory
get_target_property(ferror_LibLocation ferror LOCATION)
add_custom_command(TARGET quasi_newton_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${ferror_LibLocation} $<TARGET_FILE_DIR:quasi_newton_example>
)

get_target_property(linalg_LibLocation linalg LOCATION)
add_custom_command(TARGET quasi_newton_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${linalg_LibLocation} $<TARGET_FILE_DIR:quasi_newton_example>
)

get_target_property(blas_LibLocation blas LOCATION)
add_custom_command(TARGET quasi_newton_example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${blas_LibLocation} $<TARGET_FILE_DIR:quasi_newton_example>
)
