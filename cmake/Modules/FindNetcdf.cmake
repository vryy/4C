#
# Find the NETCDF includes and libraries
#
# NETCDF_INCLUDE_DIR - where to find autopack.h
# NETCDF_LIBRARIES   - List of fully qualified libraries to link against.
# NETCDF_FOUND       - Do not attempt to use if "no" or undefined.

IF (NETCDF_INCLUDE_DIR)
  # Already in cache, be silent
  SET(NETCDF_FIND_QUIETLY TRUE)
ENDIF (NETCDF_INCLUDE_DIR)

FIND_PATH(NETCDF_INCLUDE_DIR netcdf.h
  ${INCLUDE_INSTALL_DIR}
  /usr/include/netcdf
  /usr/include/netcdf-3
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(NETCDF_LIBRARY NAMES libnetcdf.a libnetcdf.so
  HINTS
  ${LIB_INSTALL_DIR}
  /usr/lib/netcdf-3
  /usr/local/lib
  /usr/lib
)

MARK_AS_ADVANCED(NETCDF_INCLUDE_DIR NETCDF_LIBRARY)

# Per-recommendation
SET( NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
SET( NETCDF_LIBRARIES    ${NETCDF_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NETCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDE_DIRS)
