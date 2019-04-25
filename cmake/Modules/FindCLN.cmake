#
# Find the CLN includes and libraries
#
# CLN_INCLUDE_DIR - where to find autopack.h
# CLN_LIBRARIES   - List of fully qualified libraries to link against.
# CLN_FOUND       - Do not attempt to use if "no" or undefined.

IF (CLN_INCLUDE_DIR)
  # Already in cache, be silent
  SET(CLN_FIND_QUIETLY TRUE)
ENDIF (CLN_INCLUDE_DIR)

FIND_PATH(CLN_INCLUDE_DIR cln/cln.h
  ${INCLUDE_INSTALL_DIR}
  ${SEARCH_CLN_INCLUDE_DIR}
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(CLN_LIBRARY libcln.so
  ${LIB_INSTALL_DIR}
  ${SEARCH_CLN_LIBRARY_DIR}
  /usr/local/lib
  /usr/lib
  /usr/lib64
)

MESSAGE("II use CLN in ${CLN_LIBRARY}")

SET( CLN_INCLUDE_DIRS ${CLN_INCLUDE_DIR})
SET( CLN_LIBRARIES    ${CLN_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CLN DEFAULT_MSG CLN_LIBRARIES CLN_INCLUDE_DIRS)

MARK_AS_ADVANCED(CLN_INCLUDE_DIRS CLN_LIBRARIES)