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
