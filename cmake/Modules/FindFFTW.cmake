#
# Find the FFTW includes and libraries
#
# FFTW_INCLUDE_DIR - where to find autopack.h
# FFTW_LIBRARIES   - List of fully qualified libraries to link against.
# FFTW_FOUND       - Do not attempt to use if "no" or undefined.

if(FFTW_INCLUDE_DIR)
  # Already in cache, be silent
  set(FFTW_FIND_QUIETLY TRUE)
endif(FFTW_INCLUDE_DIR)

find_path(
  FFTW_INCLUDE_DIR
  fftw3.h
  ${INCLUDE_INSTALL_DIR}
  ${FFTW_INCLUDE_DIR}
  /usr/local/include
  /usr/include
  )

find_library(
  FFTW_LIBRARY
  fftw3
  ${LIB_INSTALL_DIR}
  ${FFTW_LIBRARY_DIR}
  /usr/local/lib
  /usr/lib
  /usr/lib64
  )

set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(FFTW_LIBRARIES ${FFTW_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_DIRS)

mark_as_advanced(FFTW_INCLUDE_DIRS FFTW_LIBRARIES)
