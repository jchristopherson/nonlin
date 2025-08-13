# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    nonlin
    LANGUAGES Fortran
    VERSION 2.0.1
)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED TRUE)

# Confgiure everything
add_subdirectory(configure)

# Source
add_subdirectory(src)

# Installation Instructions
install(
  EXPORT ${PROJECT_NAME}-targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

# Examples
option(BUILD_NONLIN_EXAMPLES "Build nonlin examples")
message(STATUS "Build nonlin examples: ${BUILD_NONLIN_EXAMPLES}")
if (BUILD_NONLIN_EXAMPLES)
  add_subdirectory(examples)
endif()