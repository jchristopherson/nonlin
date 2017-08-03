# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.0)
project(nonlin C CXX Fortran)

# Define version information
set(nonlin_VERSION_MAJOR 1)
set(nonlin_VERSION_MINOR 0)
set(nonlin_VERSION_PATCH 0)

# Define output directories, if undefined
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  message(STATUS "Library output directories undefined.  Using default directories.")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    # Debug Build
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug/mod)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug/mod)
    endif ()
  elseif (CMAKE_BUILD_TYPE MATCHES Release)
    # Release Build
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release/mod)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release/mod)
    endif ()
  else ()
    # Default Condition
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug/mod)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug/mod)
    endif ()
  endif ()
endif ()

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  # Options: https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Specify source files
file(GLOB_RECURSE sources src/*.f90)

# Define necessary CMake variables, and define the dependency.
if (BUILD_SHARED_LIBS)
  set(ExternalProjectCMakeArgs
      -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
      -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
      -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
      -CMAKE_Fortran_MODULE_DIRECTORY=${CMAKE_Fortran_MODULE_DIRECTORY}
      -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
  add_subdirectory(external/linalg)
  # add_subdirectory(external/linalg/ferror)
else()
  # Building a static library, simply add source code to this library
  file(GLOB_RECURSE linalg_sources 
    external/linalg/src/*.f90 
    external/linalg/LAPACK/BLAS/src/*.f
    external/linalg/LAPACK/src/*.f 
    external/linalg/src/QRUPDATE/*.f
    external/linalg/ferror/src/*.f90)
  list(APPEND sources ${linalg_sources})
endif()

# If BUILD_SHARED_LIBS is defined, then a shared library (DLL) will be built;
# else, a static library will be built.
add_library(nonlin ${sources})
if (BUILD_SHARED_LIBS)
  target_link_libraries(nonlin linalg ferror)
endif ()

# ------------------------------------------------------------------------------
# EXAMPLES
# ------------------------------------------------------------------------------
option(BUILD_EXAMPLES "Build examples?" OFF)
if (BUILD_EXAMPLES)
  # Inform the user we're building the examples
  message(STATUS "Building NONLIN examples.")

  # Quasi-Newton Example
  add_executable(quasi_newton_example examples/nonlin_quasi_newton_example.f90)
  target_link_libraries(quasi_newton_example nonlin)

  # Least-Squares Example
  add_executable(least_squares_example 
    examples/nonlin_least_squares_example.f90)
  target_link_libraries(least_squares_example nonlin)

  # Polynomial Fitting Example
  add_executable(poly_fit_example examples/nonlin_poly_fit_example.f90)
  target_link_libraries(poly_fit_example nonlin)

  # Nelder-Mead Example
  add_executable(nm_opt_example examples/nonlin_nm_opt_example.f90)
  target_link_libraries(nm_opt_example nonlin)

# --------------------
endif()

# ------------------------------------------------------------------------------
# TESTS
# ------------------------------------------------------------------------------
set(nonlin_test_sources external/linalg/tests/test_core.f90
  tests/nonlin_test.f90 tests/nonlin_test_jacobian.f90
  tests/nonlin_test_solve.f90 tests/nonlin_test_poly.f90
  tests/nonlin_test_optimize.f90)
add_executable(nonlin_test ${nonlin_test_sources})
target_link_libraries(nonlin_test nonlin)
add_custom_command(OUTPUT nonlin_tests DEPENDS ${nonlin_test_sources}
  COMMAND nonlin_test)
add_custom_target(run_nonlin_tests ALL DEPENDS nonlin_tests)

# C Test
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests
  ${PROJECT_SOURCE_DIR}/external/linalg/ferror/include
  ${PROJECT_SOURCE_DIR}/external/linalg/tests 
  ${PROJECT_SOURCE_DIR}/external/linalg/include)
set(c_nonlin_test_sources tests/nonlin_c_test.c tests/nonlin_c_test_solve.c
  external/linalg/tests/c_test_core.c)
add_executable(nonlin_c_test ${c_nonlin_test_sources})
target_link_libraries(nonlin_c_test nonlin m)
add_custom_command(OUTPUT nonlin_c_tests DEPENDS ${c_nonlin_test_sources}
  COMMAND nonlin_c_test)
add_custom_target(run_nonlin_c_tests ALL DEPENDS nonlin_c_tests)
