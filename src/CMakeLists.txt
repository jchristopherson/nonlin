# Define the NONLIN source files
set(nonlin_sources
    nonlin_least_squares.f90
    nonlin_linesearch.f90
    nonlin_optimize.f90
    nonlin_polynomials.f90
    nonlin_solve.f90
    nonlin_constants.f90
    nonlin_core.f90
    nonlin_vecfcn_helper.f90
    nonlin_fcn1var_helper.f90
    nonlin_fcnnvar_helper.f90
    nonlin_equation_solver.f90
    nonlin_equation_solver_1var.f90
    nonlin_equation_optimizer.f90
    nonlin_optimize_nelder_mead.f90
    nonlin_optimize_line_search.f90
    nonlin_optimize_bfgs.f90
    nonlin_solve_line_search.f90
    nonlin_solve_quasi_newton.f90
    nonlin_solve_newton.f90
    nonlin_solve_brent.f90
    nonlin_solve_newton1var.f90
)

# Build the library
add_library(nonlin ${nonlin_sources})
target_link_libraries(nonlin ${LINALG_LIBRARIES} ${FERROR_LIBRARIES})

# ------------------------------------------------------------------------------
# INSTALLATION INSTRUCTIONS
# ------------------------------------------------------------------------------
# Define target information
set_property(TARGET nonlin PROPERTY VERSION ${NONLIN_VERSION})
set_property(TARGET nonlin PROPERTY SOVERSION ${NONLIN_MAJOR_VERSION})
set_property(TARGET nonlin PROPERTY INTERFACE_nonlin_MAJOR_VERSION ${NONLIN_MAJOR_VERSION})
set_property(TARGET nonlin APPEND PROPERTY COMPATIBLE_INTERFACE_STRING nonlin_MAJOR_VERSION)

# Locate the "include" directories
set(nonlin_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

# Define the installation instructions
install(TARGETS nonlin EXPORT nonlinTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    CONFIGURATIONS Release
    INCLUDES DESTINATION include
)
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    COMPONENT Devel
)

# Include the documentation
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION doc)
install(FILES ${PROJECT_SOURCE_DIR}/doc/refman.pdf DESTINATION doc)

# Define the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/nonlinConfigVersion.cmake"
    VERSION ${NONLIN_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(EXPORT nonlinTargets
    FILE "${CMAKE_BINARY_DIR}/nonlinTargets.cmake"
)

# Define the configuration file
configure_file(
    "${PROJECT_SOURCE_DIR}/nonlinConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/nonlinConfig.cmake"
    COPYONLY
)

set(ConfigPackageLocation lib/cmake/nonlin)
install(
    EXPORT nonlinTargets
    FILE nonlinTargets.cmake
    DESTINATION ${ConfigPackageLocation}
)
install(
    FILES
        "${CMAKE_BINARY_DIR}/nonlinConfig.cmake"
        "${CMAKE_BINARY_DIR}/nonlinConfigVersion.cmake"
    DESTINATION
        ${ConfigPackageLocation}
    COMPONENT
        Devel
)
