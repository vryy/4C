#
# Find the TETGEN includes and libraries
#
# TETGEN_INCLUDE_DIR - where to find autopack.h
# TETGEN_LIBRARIES   - List of fully qualified libraries to link against.
# TETGEN_FOUND       - Do not attempt to use if "no" or undefined.

IF (TETGEN_INCLUDE_DIR)
  # Already in cache, be silent
  SET(TETGEN_FIND_QUIETLY TRUE)
ENDIF (TETGEN_INCLUDE_DIR)

FIND_PATH(TETGEN_INCLUDE_DIR tetgen.h
  ${INCLUDE_INSTALL_DIR}
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(TETGEN_LIBRARY NAMES tetgen
  PATHS
  ${LIB_INSTALL_DIR}
  /usr/local/lib
  /usr/lib
)

SET( TETGEN_INCLUDE_DIRS ${TETGEN_INCLUDE_DIR})
SET( TETGEN_LIBRARIES    ${TETGEN_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TETGEN DEFAULT_MSG TETGEN_LIBRARIES TETGEN_INCLUDE_DIRS)

MARK_AS_ADVANCED(TETGEN_INCLUDE_DIRS TETGEN_LIBRARIES)
