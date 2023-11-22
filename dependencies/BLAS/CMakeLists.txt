# If found, use the installed version; else, import the library
if (${BLAS_FOUND})
    # Inform the user of what's going on
    message(STATUS "BLAS library found.")
else()
    message(STATUS "BLAS library not found.  The reference BLAS will be used.")
endif()