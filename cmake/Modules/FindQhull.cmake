#
# Find the QHULL includes and libraries
#
# QHULL_INCLUDE_DIR - where to find autopack.h
# QHULL_LIBRARIES   - List of fully qualified libraries to link against.
# QHULL_FOUND       - Do not attempt to use if "no" or undefined.

if(QHULL_INCLUDE_DIR)
  # Already in cache, be silent
  set(QHULL_FIND_QUIETLY TRUE)
endif(QHULL_INCLUDE_DIR)

find_path(
  QHULL_INCLUDE_DIR
  qhull/qhull.h
  ${INCLUDE_INSTALL_DIR}
  ${QHULL_INCLUDE_DIR}
  /usr/local/include
  /usr/include
  )

find_library(
  QHULL_LIBRARY
  libqhull.so
  ${LIB_INSTALL_DIR}
  ${QHULL_LIBRARY_DIR}
  NO_DEFAULT_PATH
  )

message("II use QHULL in ${QHULL_LIBRARY}")

set(QHULL_INCLUDE_DIRS ${QHULL_INCLUDE_DIR})
set(QHULL_LIBRARIES ${QHULL_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qhull DEFAULT_MSG QHULL_LIBRARIES QHULL_INCLUDE_DIRS)

mark_as_advanced(QHULL_INCLUDE_DIRS QHULL_LIBRARIES)
