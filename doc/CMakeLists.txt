#
# Build routine/test for doxygen documentation
option(BUILD_DOXYGEN "Build doxygen documentation" ON)
if(BUILD_DOXYGEN)

  # add doc-specific cmake modules
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

  find_package(Doxygen REQUIRED)

  # set input and output files
  set(DOXYGEN_CONFIG_IN ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in)
  set(DOXYGEN_OUT_DIRECTORY ${CMAKE_BINARY_DIR}/doc/)
  set(DOXYGEN_CONFIG_OUT ${DOXYGEN_OUT_DIRECTORY}/Doxyfile)

  # find mathjax
  message(STATUS "Doxygen: Use mathjax for math equations")
  find_package(MathJax REQUIRED)

  # configure doxygen configuration
  configure_file(${DOXYGEN_CONFIG_IN} ${DOXYGEN_CONFIG_OUT} @ONLY)

  add_custom_target(
    documentation ALL
    COMMAND Doxygen::doxygen ${DOXYGEN_CONFIG_OUT}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    VERBATIM
    COMMENT "Build doxygen"
    )

  add_dependencies(documentation setup_mathjax)

  list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
endif()
