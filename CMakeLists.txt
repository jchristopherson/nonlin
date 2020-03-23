# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.7)
project(nonlin Fortran)

# Define version information
set(NONLIN_MAJOR_VERSION 1)
set(NONLIN_MINOR_VERSION 4)
set(NONLIN_PATCH_VERSION 0)
set(NONLIN_VERSION ${NONLIN_MAJOR_VERSION}.${NONLIN_MINOR_VERSION}.${NONLIN_PATCH_VERSION})

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# By default, shared library
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Get compiler info
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# Export all symbols on Windows when building shared libraries
SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Locate the module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/include)

# Define output directories, if undefined
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    message(STATUS "NONLIN library output directories undefined.  Using default directories.")
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        # Debug Build
        if (BUILD_SHARED_LIBS)
            set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
            set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
        else()
            set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
            set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
        endif()
    elseif (CMAKE_BUILD_TYPE MATCHES Release)
        # Release Build
        if (BUILD_SHARED_LIBS)
            set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
            set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
        else()
            set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
            set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
        endif()
    else()
        # Default Condition
        if (BUILD_SHARED_LIBS)
            set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
            set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
        else()
            set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
            set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
        endif()
    endif()
endif()

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -Wl,--allow-multiple-definition")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall -Wno-c-binding-type -Wl,--allow-multiple-definition")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -Wl,--allow-multiple-definition")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall -Wl,--allow-multiple-definition")
endif ()

# Locate Dependencies
find_package(ferror 1.3.0)
find_package(linalg 1.6.0)

if (ferror_FOUND)
  message(STATUS "FERROR library found.")
  set(ferror_LIBRARIES ferror)
endif()

if (linalg_FOUND)
    message(STATUS "LINALG library found.")
    set(linalg_LIBRARIES linalg)
endif()

# Locate the source directory
add_subdirectory(src)

# ------------------------------------------------------------------------------
# EXAMPLES
# ------------------------------------------------------------------------------
option(BUILD_NONLIN_EXAMPLES "Build NONLIN examples?" OFF)
if (BUILD_NONLIN_EXAMPLES)
  # Inform the user we're building the examples
  message(STATUS "Building NONLIN examples.")

  # Build the examples
  add_subdirectory(examples)
endif()

# ------------------------------------------------------------------------------
# TESTS
# ------------------------------------------------------------------------------
option(BUILD_NONLIN_TESTS "Build NONLIN tests?" OFF)
if (BUILD_NONLIN_TESTS)
    # Inform the user we're building the tests
    message(STATUS "Building NONLIN tests.")

    # Build the tests
    add_subdirectory(tests)
endif()
