find_package(FFTW REQUIRED)

if(FFTW_FOUND)
  message(STATUS "FFTW include directory: ${FFTW_INCLUDE_DIR}")
  message(STATUS "FFTW libraries: ${FFTW_LIBRARY}")

  target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE fftw::fftw)
endif()
