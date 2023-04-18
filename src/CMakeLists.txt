# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the NONLIN source files
set(NONLIN_SOURCES
    ${dir}/nonlin_least_squares.f90
    ${dir}/nonlin_linesearch.f90
    ${dir}/nonlin_optimize.f90
    ${dir}/nonlin_polynomials.f90
    ${dir}/nonlin_solve.f90
    ${dir}/nonlin_constants.f90
    ${dir}/nonlin_core.f90
    ${dir}/nonlin_vecfcn_helper.f90
    ${dir}/nonlin_fcn1var_helper.f90
    ${dir}/nonlin_fcnnvar_helper.f90
    ${dir}/nonlin_equation_solver.f90
    ${dir}/nonlin_equation_solver_1var.f90
    ${dir}/nonlin_equation_optimizer.f90
    ${dir}/nonlin_optimize_nelder_mead.f90
    ${dir}/nonlin_optimize_line_search.f90
    ${dir}/nonlin_optimize_bfgs.f90
    ${dir}/nonlin_solve_line_search.f90
    ${dir}/nonlin_solve_quasi_newton.f90
    ${dir}/nonlin_solve_newton.f90
    ${dir}/nonlin_solve_brent.f90
    ${dir}/nonlin_solve_newton1var.f90
)
set(NONLIN_SOURCES ${NONLIN_SOURCES} PARENT_SCOPE)