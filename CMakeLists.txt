# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    nonlin
    LANGUAGES Fortran
    VERSION 2.0.0
)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
find_package(BLAS)
find_package(LAPACK)
add_subdirectory(dependencies)

if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  message(STATUS "BLAS/LAPACK could not be found.  A reference version will be employed.")
  include(FetchContent)
  FetchContent_Declare(
      lapack
      GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
  )
  FetchContent_MakeAvailable(lapack)
  set(BLAS_LIBRARIES blas)
  set(LAPACK_LIBRARIES lapack)
endif()

# Source
add_subdirectory(src)
add_fortran_library(
    ${PROJECT_NAME}
    ${PROJECT_INCLUDE_DIR}
    ${CMAKE_INSTALL_INCLUDEDIR}
    ${PROJECT_VERSION}
    ${PROJECT_VERSION_MAJOR}
    ${NONLIN_SOURCES}
)
target_link_libraries(
  ${PROJECT_NAME} 
  ${BLAS_LIBRARIES} 
  ${LAPACK_LIBRARIES}
)
link_library(${PROJECT_NAME} ${ferror_LIBRARY} ${ferror_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${linalg_LIBRARY} ${linalg_INCLUDE_DIR})

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