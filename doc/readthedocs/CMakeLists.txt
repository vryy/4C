#
# Build routine/test for readthedocs documentation
option(BACI_BUILD_READTHEDOCS "Build readthedocs documentation" OFF)
if(BACI_BUILD_READTHEDOCS)

  # add doc-specific cmake modules
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

  # set input and output files and directories
  set(BACI_READTHEDOCS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/doc/readthedocs)
  set(BACI_READTHEDOCS_CONFIG_DIR ${CMAKE_SOURCE_DIR}/doc/readthedocs-config)
  set(BACI_READTHEDOCS_VENV_DIR ${BACI_READTHEDOCS_CONFIG_DIR}/baci-python-rtd-venv)
  set(BACI_READTHEDOCS_CONFIG_IN ${BACI_READTHEDOCS_CONFIG_DIR}/conf.py.in)
  set(BACI_READTHEDOCS_OUT_DIR ${CMAKE_BINARY_DIR}/doc/readthedocs)
  set(BACI_READTHEDOCS_CONFIG_OUT ${BACI_READTHEDOCS_OUT_DIR}/conf.py)
  set(BACI_READTHEDOCS_HTML_DIR ${BACI_READTHEDOCS_OUT_DIR}/html)

  file(MAKE_DIRECTORY ${BACI_READTHEDOCS_OUT_DIR}/_static)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/tmp)
  file(REMOVE_RECURSE "${BACI_READTHEDOCS_VENV_DIR}")

  #install the python virtual environment with sphinx
  find_package(Sphinx REQUIRED)

  if(NOT DEFINED SPHINX_THEME)
    set(SPHINX_THEME sphinx_rtd_theme)
  endif()

  # configure readthedocs configuration
  configure_file(${BACI_READTHEDOCS_CONFIG_IN} ${BACI_READTHEDOCS_CONFIG_OUT} @ONLY)

  add_custom_target(
    readthedocs ALL
    COMMAND create_rtd
    COMMAND ${baciname} -h > baci-help.txt
    COMMAND
      ${SPHINX_EXECUTABLE} -q -b html # -W turns warnings into errors!
      -i "${CMAKE_CURRENT_SOURCE_DIR}" -i "${BACI_READTHEDOCS_OUT_DIR}" -i
      "${CMAKE_SOURCE_DIR}/tests/framework-test" -o "${BACI_READTHEDOCS_HTML_DIR}" -c
      "${BACI_READTHEDOCS_OUT_DIR}" -s "${CMAKE_BINARY_DIR}/doc/tmp"
    DEPENDS create_rtd ${baciname}
    COMMENT "Build readthedocs documentation with Sphinx"
    )

  list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
endif()
