#
# Build routine/test for readthedocs documentation
option(BACI_BUILD_READTHEDOCS "Build readthedocs documentation" OFF)
if(BACI_BUILD_READTHEDOCS)

  # add doc-specific cmake modules
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

  # set input and output files and directories
  set(BACI_VENV_DIR ${CMAKE_SOURCE_DIR}/utilities/baci-python-venv)
  set(READTHEDOCS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/doc/readthedocs)
  set(READTHEDOCS_CONFIG_DIR ${CMAKE_SOURCE_DIR}/doc/readthedocs-config)
  set(READTHEDOCS_CONFIG_IN ${READTHEDOCS_CONFIG_DIR}/conf.py.in)
  set(READTHEDOCS_OUT_DIR ${CMAKE_BINARY_DIR}/doc/readthedocs)
  set(READTHEDOCS_CONFIG_OUT ${READTHEDOCS_OUT_DIR}/conf.py)
  set(READTHEDOCS_HTML_DIR ${READTHEDOCS_OUT_DIR}/html)

  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/tmp)

  #install the python virtual environment with sphinx
  find_package(Sphinx REQUIRED)

  if(NOT DEFINED BACI_SPHINX_THEME)
    set(BACI_SPHINX_THEME sphinx_rtd_theme)
  endif()

  # configure readthedocs configuration
  configure_file(${READTHEDOCS_CONFIG_IN} ${READTHEDOCS_CONFIG_OUT} @ONLY)
  message("Sphinx configure file conf.in is written based on ${READTHEDOCS_CONFIG_IN}")

  add_custom_target(
    readthedocs ALL
    COMMAND create_rtd
    COMMAND ${baciname} -h > baci-help.txt
    COMMAND
      ${SPHINX_EXECUTABLE} -q -b html # -W turns warnings into errors!
      -i "${CMAKE_CURRENT_SOURCE_DIR}" -i "${READTHEDOCS_OUT_DIR}" -i
      "${CMAKE_SOURCE_DIR}/tests/framework-test" -o "${READTHEDOCS_HTML_DIR}" -c
      "${READTHEDOCS_OUT_DIR}" -s "${CMAKE_BINARY_DIR}/doc/tmp"
    DEPENDS create_rtd ${baciname}
    COMMENT "Build readthedocs documentation with Sphinx"
    )

  list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
endif()
