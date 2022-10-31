#
# Build routine/test for doxygen documentation
option(BUILD_DOXYGEN "Build doxygen documentation" ON)

# prefer doxygen at specified path instead of default path
find_program(
  DOXYGEN doxygen
  PATHS ${DOXYGEN_BASE_PATH}
  NO_DEFAULT_PATH
  )
find_program(DOXYGEN doxygen)
if(DOXYGEN)
  message("Found doxygen: ${DOXYGEN}")

  # set input and output files
  set(DOXYGEN_CONFIG_IN ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in)
  set(DOXYGEN_OUT_DIRECTORY ${CMAKE_BINARY_DIR}/doc/)
  set(DOXYGEN_CONFIG_OUT ${DOXYGEN_OUT_DIRECTORY}/Doxyfile)

  # find mathjax
  include(FindMathJax.cmake)

  # configure doxygen configuration
  configure_file(${DOXYGEN_CONFIG_IN} ${DOXYGEN_CONFIG_OUT} @ONLY)

  add_custom_target(
    documentation ALL
    COMMAND ${DOXYGEN} ${DOXYGEN_CONFIG_OUT}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    VERBATIM
    )

  if(DOXYGEN_USE_LOCAL_MATHJAX)
    add_custom_command(
      TARGET documentation
      POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_OUT_DIRECTORY}/html/${DOXYGEN_MATHJAX_RELPATH}
      COMMAND
        ${CMAKE_COMMAND} -E copy_directory ${DOXYGEN_LOCAL_MATHJAX_BASEPATH}
        ${DOXYGEN_OUT_DIRECTORY}/html/${DOXYGEN_MATHJAX_RELPATH}
      )
  endif()
else()
  message("Doxygen not found.")
endif()
