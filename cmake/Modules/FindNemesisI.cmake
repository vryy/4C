#
# Find the NEMESIS includes and libraries
#
# NEMESIS_INCLUDE_DIR - where to find autopack.h
# NEMESIS_LIBRARIES   - List of fully qualified libraries to link against.
# NEMESIS_FOUND       - Do not attempt to use if "no" or undefined.

IF (NEMESIS_INCLUDE_DIR)
  # Already in cache, be silent
  SET(NEMESIS_FIND_QUIETLY TRUE)
ENDIF (NEMESIS_INCLUDE_DIR)

FIND_PATH(NEMESIS_INCLUDE_DIR ne_nemesisI.h
  ${INCLUDE_INSTALL_DIR}
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(NEMESIS_LIBRARY libnemIc.a
  ${LIB_INSTALL_DIR}
  /usr/local/lib
  /usr/lib
)

SET( NEMESIS_INCLUDE_DIRS ${NEMESIS_INCLUDE_DIR})
SET( NEMESIS_LIBRARIES    ${NEMESIS_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NEMESIS DEFAULT_MSG NEMESIS_LIBRARIES NEMESIS_INCLUDE_DIRS)

MARK_AS_ADVANCED(NEMESIS_INCLUDE_DIRS NEMESIS_LIBRARIES)
