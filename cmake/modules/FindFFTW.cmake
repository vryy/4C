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

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARY FFTW_INCLUDE_DIR)

if(FFTW_FOUND AND NOT TARGET fftw::fftw)
  add_library(fftw::fftw UNKNOWN IMPORTED)
  set_target_properties(
    fftw::fftw
    PROPERTIES IMPORTED_LOCATION "${FFTW_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
    )
endif()

if(FFTW_FOUND)
  message(STATUS "FFTW include directory: ${FFTW_INCLUDE_DIR}")
  message(STATUS "FFTW libraries: ${FFTW_LIBRARY}")

  list(APPEND BACI_ALL_ENABLED_EXTERNAL_LIBS fftw::fftw)
  target_compile_definitions(fftw::fftw INTERFACE "HAVE_FFTW")
  set(HAVE_FFTW ON)
else()
  message(WARNING "FFTW not found. no FFTW support in BACI.")
  message(WARNING "You cannot run simulations that need a discrete FFT.")
endif()
