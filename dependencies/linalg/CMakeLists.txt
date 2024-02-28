include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    linalg
    GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
)

FetchContent_MakeAvailable(linalg)

set(linalg_INCLUDE_DIR ${linalg_BINARY_DIR}/include)
set(linalg_INCLUDE_DIR ${linalg_INCLUDE_DIR} PARENT_SCOPE)
configure_file(
    "${linalg_SOURCE_DIR}/include/linalg.h"
    "${linalg_INCLUDE_DIR}/linalg.h"
    COPYONLY
)
set(linalg_LIBRARY linalg)
set(linalg_LIBRARY ${linalg_LIBRARY} PARENT_SCOPE)
