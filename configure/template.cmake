@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(ferror QUIET)
find_dependency(linalg QUIET)


if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()