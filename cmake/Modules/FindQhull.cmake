
#
# Find the QHULL includes and libraries
#
# QHULL_INCLUDE_DIR - where to find autopack.h
# QHULL_LIBRARIES   - List of fully qualified libraries to link against.
# QHULL_FOUND       - Do not attempt to use if "no" or undefined.

IF (QHULL_INCLUDE_DIR)
  # Already in cache, be silent
  SET(QHULL_FIND_QUIETLY TRUE)
ENDIF (QHULL_INCLUDE_DIR)

FIND_PATH(QHULL_INCLUDE_DIR qhull/qhull.h
  ${INCLUDE_INSTALL_DIR}
  ${QHULL_INCLUDE_DIR}
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(QHULL_LIBRARY libqhull.so
  ${LIB_INSTALL_DIR}
  ${QHULL_LIBRARY_DIR}
  NO_DEFAULT_PATH
)

MESSAGE("II use QHULL in ${QHULL_LIBRARY}")

SET( QHULL_INCLUDE_DIRS ${QHULL_INCLUDE_DIR})
SET( QHULL_LIBRARIES    ${QHULL_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Qhull DEFAULT_MSG QHULL_LIBRARIES QHULL_INCLUDE_DIRS)

MARK_AS_ADVANCED(QHULL_INCLUDE_DIRS QHULL_LIBRARIES)
