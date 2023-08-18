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

message(STATUS "Use CLN in ${CLN_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLN DEFAULT_MSG CLN_LIBRARY CLN_INCLUDE_DIR)

if(CLN_FOUND AND NOT TARGET cln::cln)
  add_library(cln::cln UNKNOWN IMPORTED)
  set_target_properties(
    cln::cln
    PROPERTIES IMPORTED_LOCATION "${CLN_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${CLN_INCLUDE_DIR}"
    )
endif()

if(CLN_FOUND)
  list(APPEND BACI_ALL_ENABLED_EXTERNAL_LIBS cln::cln)
  message(STATUS "CLN include directory: ${CLN_INCLUDE_DIR}")
  message(STATUS "CLN library directory: ${CLN_LIBRARY}")
endif()
