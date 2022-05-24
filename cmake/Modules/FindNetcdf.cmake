#
# Find the NETCDF includes and libraries
#
# NETCDF_INCLUDE_DIR - where to find autopack.h
# NETCDF_LIBRARIES   - List of fully qualified libraries to link against.
# NETCDF_FOUND       - Do not attempt to use if "no" or undefined.

if(NETCDF_INCLUDE_DIR)
  # Already in cache, be silent
  set(NETCDF_FIND_QUIETLY TRUE)
endif(NETCDF_INCLUDE_DIR)

find_path(
  NETCDF_INCLUDE_DIR
  netcdf.h
  ${INCLUDE_INSTALL_DIR}
  /usr/include/netcdf
  /usr/include/netcdf-3
  /usr/local/include
  /usr/include
  )

find_library(
  NETCDF_LIBRARY
  NAMES libnetcdf.a libnetcdf.so
  HINTS ${LIB_INSTALL_DIR} /usr/lib/netcdf-3 /usr/local/lib /usr/lib
  )

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_LIBRARY)

# Per-recommendation
set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
set(NETCDF_LIBRARIES ${NETCDF_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Netcdf DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDE_DIRS)
