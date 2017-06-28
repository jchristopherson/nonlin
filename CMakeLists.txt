# Master CMAKE Build Script
cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
project(nonlin)

# TO DO:
# Create link to linalg project, and then ensure CMAKE builds the linalg library
# for use in this project - force use of a shared library?
# Regardless, here's some info on how to accomplish:
# http://foonathan.net/blog/2016/07/07/cmake-dependency-handling.html
# https://softwareengineering.stackexchange.com/questions/114719/organizing-git-repositories-with-common-nested-sub-modules
# http://diracprogram.org/doc/release-14/programmers/external_projects.html - looks to be a good reference

# Define version information
set(nonlin_VERSION_MAJOR 1)
set(nonlin_VERSION_MINOR 0)
set(nonlin_VERSION_PATCH 0)
set(nonlin_VERSION_TWEAK 0)

# Specify output directories
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
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
  add_external(linalg "")
  add_dependencies(nonlin linalg)

  # TO DO: 
  # - Still link with ferror
  # - Be sure to check if the libraries already exist.  If so, simply link
  # - Figure out how to deal with .mod files
else()
  # Building a static library, simply add source code to this library
  file(GLOB_RECURSE linalg_sources 
    external/linalg/src/*.f90 
    external/linalg/src/BLAS/*.f
    external/linalg/src/LAPACK/*.f 
    external/linalg/src/QRUPDATE/*.f)
  list(APPEND sources ${linalg_sources})
  
  # Also include the stuff from FERROR
  file(GLOB_RECURSE ferror_sources external/ferror/src/*.f90)
  list(APPEND sources ${ferror_sources})
endif()

# Define testing source files
set(testing_sources ${sources})

# Append any files in the tests directory
list(APPEND testing_sources external/linalg/tests/test_core.f90
  tests/nonlin_test.f90 tests/nonlin_test_jacobian.f90
  tests/nonlin_test_solve.f90)

# If BUILD_SHARED_LIBS is defined, then a shared library (DLL) will be built;
# else, a static library will be built.
add_library(nonlin ${sources})

# Build the test executable
 add_executable(nonlin_test ${testing_sources})