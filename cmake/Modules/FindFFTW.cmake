#
# Find the FFTW includes and libraries
#
# FFTW_INCLUDE_DIR - where to find autopack.h
# FFTW_LIBRARIES   - List of fully qualified libraries to link against.
# FFTW_FOUND       - Do not attempt to use if "no" or undefined.

IF (FFTW_INCLUDE_DIR)
  # Already in cache, be silent
  SET(FFTW_FIND_QUIETLY TRUE)
ENDIF (FFTW_INCLUDE_DIR)

FIND_PATH(FFTW_INCLUDE_DIR fftw3.h
  ${INCLUDE_INSTALL_DIR}
  ${FFTW_INCLUDE_DIR}
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(FFTW_LIBRARY fftw3
  ${LIB_INSTALL_DIR}
  ${FFTW_LIBRARY_DIR}
  /usr/local/lib
  /usr/lib
  /usr/lib64
)

SET( FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
SET( FFTW_LIBRARIES    ${FFTW_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_DIRS)

MARK_AS_ADVANCED(FFTW_INCLUDE_DIRS FFTW_LIBRARIES)
