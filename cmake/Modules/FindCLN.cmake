#
# Find the CLN includes and libraries
#
# CLN_INCLUDE_DIR - where to find autopack.h
# CLN_LIBRARIES   - List of fully qualified libraries to link against.
# CLN_FOUND       - Do not attempt to use if "no" or undefined.

if(CLN_INCLUDE_DIR)
  # Already in cache, be silent
  set(CLN_FIND_QUIETLY TRUE)
endif(CLN_INCLUDE_DIR)

find_path(
  CLN_INCLUDE_DIR
  cln/cln.h
  ${INCLUDE_INSTALL_DIR}
  ${SEARCH_CLN_INCLUDE_DIR}
  /usr/local/include
  /usr/include
  )

find_library(
  CLN_LIBRARY
  libcln.so
  ${LIB_INSTALL_DIR}
  ${SEARCH_CLN_LIBRARY_DIR}
  /usr/local/lib
  /usr/lib
  /usr/lib64
  )

message("II use CLN in ${CLN_LIBRARY}")

set(CLN_INCLUDE_DIRS ${CLN_INCLUDE_DIR})
set(CLN_LIBRARIES ${CLN_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLN DEFAULT_MSG CLN_LIBRARIES CLN_INCLUDE_DIRS)

mark_as_advanced(CLN_INCLUDE_DIRS CLN_LIBRARIES)
