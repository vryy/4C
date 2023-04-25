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

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qhull DEFAULT_MSG QHULL_LIBRARY QHULL_INCLUDE_DIR)

if(QHULL_FOUND AND NOT TARGET qhull::qhull)
  add_library(qhull::qhull UNKNOWN IMPORTED)
  set_target_properties(
    qhull::qhull
    PROPERTIES IMPORTED_LOCATION "${QHULL_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${QHULL_INCLUDE_DIR}"
    )
endif()

if(QHULL_FOUND)
  message(STATUS "QHULL include directory: ${QHULL_INCLUDE_DIR}")
  message(STATUS "QHULL library directory: ${QHULL_LIBRARY}")
  list(APPEND BACI_ALL_ENABLED_EXTERNAL_LIBS qhull::qhull)
endif(QHULL_FOUND)
