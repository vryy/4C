#
# Find the EXODUS includes and libraries
#
# EXODUS_INCLUDE_DIR - where to find autopack.h
# EXODUS_LIBRARIES   - List of fully qualified libraries to link against.
# EXODUS_FOUND       - Do not attempt to use if "no" or undefined.

IF (EXODUS_INCLUDE_DIR)
  # Already in cache, be silent
  SET(EXODUS_FIND_QUIETLY TRUE)
ENDIF (EXODUS_INCLUDE_DIR)

FIND_PATH(EXODUS_INCLUDE_DIR exodusII.h
  ${INCLUDE_INSTALL_DIR}
  /usr/include/netcdf
  /usr/include/netcdf-3
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(EXODUS_LIBRARY libexoIIv2c.a
  ${LIB_INSTALL_DIR}
  /usr/local/lib
  /usr/lib
)

SET( EXODUS_INCLUDE_DIRS ${EXODUS_INCLUDE_DIR})
SET( EXODUS_LIBRARIES    ${EXODUS_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(EXODUS DEFAULT_MSG EXODUS_LIBRARIES EXODUS_INCLUDE_DIRS)

MARK_AS_ADVANCED(EXODUS_INCLUDE_DIRS EXODUS_LIBRARIES)
