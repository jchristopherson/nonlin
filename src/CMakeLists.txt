# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the NONLIN source files
set(NONLIN_SOURCES
    ${dir}/nonlin.f90
    ${dir}/nonlin_error_handling.f90
    ${dir}/nonlin_multi_eqn_mult_var.f90
    ${dir}/nonlin_single_var.f90
    ${dir}/nonlin_multi_var.f90
    ${dir}/nonlin_types.f90
    ${dir}/nonlin_least_squares.f90
    ${dir}/nonlin_linesearch.f90
    ${dir}/nonlin_optimize.f90
    ${dir}/nonlin_solve.f90
    ${dir}/nonlin_helper.f90
    ${dir}/nonlin_polynomials.f90
)
set(NONLIN_SOURCES ${NONLIN_SOURCES} PARENT_SCOPE)