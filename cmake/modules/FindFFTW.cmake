# Finder for FFTW
# Exports fftw::fftw as an imported target
# Note: The FFTW_ROOT variable is automatically considered by the find_ calls below.

find_path(FFTW_INCLUDE_DIR fftw3.h)

find_library(FFTW_LIBRARY fftw3)

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

  target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE fftw::fftw)
endif()
