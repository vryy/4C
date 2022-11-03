# We do not use the bundled finder for Doxygen since we want to provide custom options

# prefer doxygen at specified path instead of default path
find_program(
  _doxygen_internal doxygen
  PATHS ${DOXYGEN_BASE_PATH}
  NO_DEFAULT_PATH
  )
# fall back to a default version if the above command did not work
find_program(_doxygen_internal doxygen REQUIRED)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Doxygen DEFAULT_MSG _doxygen_internal)

# create an imported executable consistent with the modern CMake finder
add_executable(Doxygen::doxygen IMPORTED)
set_property(TARGET Doxygen::doxygen PROPERTY IMPORTED_LOCATION ${_doxygen_internal})
