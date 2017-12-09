# Locate the dependencies
add_subdirectory(external/linalg)

# Define the NONLIN source files
set(nonlin_sources nonlin_c_binding.f90 nonlin_least_squares.f90 
    nonlin_linesearch.f90 nonlin_optimize.f90 nonlin_polynomials.f90
    nonlin_solve.f90 nonlin_types.f90)

# Build the library
add_library(nonlin ${nonlin_sources})
target_link_libraries(nonlin linalg)
