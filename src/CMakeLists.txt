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

# Dependecies
include(FetchContent)

find_package(BLAS)
find_package(LAPACK)
find_package(ferror QUIET)
find_package(linalg QUIET)

if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  message(STATUS "BLAS/LAPACK could not be found.  A reference version will be employed.")
  FetchContent_Declare(
      lapack
      GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
  )
  FetchContent_MakeAvailable(lapack)
  set(BLAS_LIBRARIES blas)
  set(LAPACK_LIBRARIES lapack)
endif()

if (NOT ferror_FOUND)
    message(STATUS "FERROR could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
    set(BUILD_FERROR TRUE)
else()
    set(FERROR_LIBRARIES ferror::ferror)
    set(BUILD_FERROR FALSE)
endif()

if (NOT linalg_FOUND)
    message(STATUS "LINALG could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
    set(BUILD_LINALG TRUE)
else()
    set(LINALG_LIBRARIES linalg::linalg)
    set(BUILD_LINALG FALSE)
endif()

add_library(${PROJECT_NAME} ${NONLIN_SOURCES})
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${FERROR_LIBRARIES}
    ${LINALG_LIBRARIES}
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set(NONLIN_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
if (NOT EXISTS "${NONLIN_MOD_DIR}")
    file(MAKE_DIRECTORY "${NONLIN_MOD_DIR}")
endif()

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${NONLIN_MOD_DIR}
)
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${NONLIN_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${NONLIN_MOD_DIR} DESTINATION ${CMAKE_INSTALL_MODULEDIR})

if (${BUILD_FERROR} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${FERROR_LIBRARIES})
endif()

if (${BUILD_LINALG} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${LINALG_LIBRARIES})
endif()
