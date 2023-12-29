#
# Build routine/test for doxygen documentation
baci_process_global_option(BACI_BUILD_DOXYGEN "Build doxygen documentation" OFF)
if(BACI_BUILD_DOXYGEN)

  # add doc-specific cmake modules
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

  find_package(Doxygen REQUIRED dot)

  # set input and output files
  set(DOXYGEN_CONFIG_IN ${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile.in)
  set(DOXYGEN_OUT_DIRECTORY ${PROJECT_BINARY_DIR}/doc/doxygen)
  set(DOXYGEN_CONFIG_OUT ${DOXYGEN_OUT_DIRECTORY}/Doxyfile)

  # find mathjax
  find_package(MathJax REQUIRED)

  # configure doxygen configuration
  configure_file(${DOXYGEN_CONFIG_IN} ${DOXYGEN_CONFIG_OUT} @ONLY)

  add_custom_target(
    doxygen ALL
    COMMAND Doxygen::doxygen ${DOXYGEN_CONFIG_OUT}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    VERBATIM
    COMMENT "Build doxygen"
    )

  add_dependencies(doxygen setup_mathjax)

  list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
endif()
