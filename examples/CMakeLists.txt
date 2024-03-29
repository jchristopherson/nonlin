# Quasi-Newton Example
add_executable(quasi_newton_example nonlin_quasi_newton_example.f90 example_problems.f90)
target_link_libraries(quasi_newton_example nonlin)

# Least-Squares Example
add_executable(least_squares_example nonlin_least_squares_example.f90 example_problems.f90)
target_link_libraries(least_squares_example nonlin)

# Polynomial Fitting Example
add_executable(poly_fit_example nonlin_poly_fit_example.f90)
target_link_libraries(poly_fit_example nonlin)

# Nelder-Mead Example
add_executable(nm_opt_example nonlin_nm_opt_example.f90 example_problems.f90)
target_link_libraries(nm_opt_example nonlin)

# Newton Example
add_executable(newton_example nonlin_newton_example.f90 example_problems.f90)
target_link_libraries(newton_example nonlin)

# Newton Example w/ Jacobian
add_executable(newton_example_jacobian nonlin_newton_solve_jacobian.f90 example_problems.f90)
target_link_libraries(newton_example_jacobian nonlin)

# Brent Example
add_executable(brent_example nonlin_brent_example.f90 example_problems.f90)
target_link_libraries(brent_example nonlin)

# BFGS Example
add_executable(bfgs_example nonlin_bfgs_example.f90 example_problems.f90)
target_link_libraries(bfgs_example nonlin)

# 1D Newton Example
add_executable(newton_1d_example nonlin_newton1d_example.f90 example_problems.f90)
target_link_libraries(newton_1d_example nonlin)

# Newton Example - 2
add_executable(newton_example_2 nonlin_newton_example_2.f90 example_problems.f90)
target_link_libraries(newton_example_2 nonlin)

# Polynomial Roots Example
add_executable(polynomial_roots_example nonlin_poly_roots_example.f90)
target_link_libraries(polynomial_roots_example nonlin)
