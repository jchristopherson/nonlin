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

# C API - 1D Newton Example
include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(c_newton_1d_example nonlin_c_newton_1d_example.c)
target_link_libraries(c_newton_1d_example nonlin)

# C API - 1D Brent Example
add_executable(c_brent_1d_example nonlin_c_brent_1d_example.c)
target_link_libraries(c_brent_1d_example nonlin)

# C API - Quasi-Newton Example
add_executable(c_quasi_newton_example nonlin_c_quasi_newton_example.c)
target_link_libraries(c_quasi_newton_example nonlin)

# C API - Newton Example
add_executable(c_newton_example nonlin_c_newton_example.c)
target_link_libraries(c_newton_example nonlin)

# C API - Least Squares Example
add_executable(c_least_squares_example nonlin_c_least_squares_example.c)
target_link_libraries(c_least_squares_example nonlin)

# C API - Nelder-Mead Example
add_executable(c_nelder_mead_example nonlin_c_nelder_mead_example.c)
target_link_libraries(c_nelder_mead_example nonlin)

# C API - BFGS Example
add_executable(c_bfgs_example nonlin_c_bfgs_example.c)
target_link_libraries(c_bfgs_example nonlin)

# C API - POLYNOMIAL Example
add_executable(c_polynomial_example nonlin_c_polynomial_example.c)
target_link_libraries(c_polynomial_example nonlin)
