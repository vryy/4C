#
# Build routine/test for readthedocs documentation
baci_process_global_option(FOUR_C_BUILD_READTHEDOCS "Build readthedocs documentation" OFF)
if(FOUR_C_BUILD_READTHEDOCS)

  # add doc-specific cmake modules
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

  # set input and output files and directories
  set(FOUR_C_VENV_DIR ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv)
  set(READTHEDOCS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/doc/readthedocs)
  set(READTHEDOCS_CONFIG_DIR ${PROJECT_SOURCE_DIR}/doc/readthedocs-config)
  set(READTHEDOCS_TPLDEPEND_DIR ${PROJECT_SOURCE_DIR}/dependencies/current)
  set(READTHEDOCS_CONFIG_IN ${READTHEDOCS_CONFIG_DIR}/conf.py.in)
  set(READTHEDOCS_OUT_DIR ${PROJECT_BINARY_DIR}/doc/readthedocs)
  set(READTHEDOCS_CONFIG_OUT ${READTHEDOCS_OUT_DIR}/conf.py)
  set(READTHEDOCS_HTML_DIR ${READTHEDOCS_OUT_DIR}/html)

  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/doc/tmp)
  file(MAKE_DIRECTORY ${READTHEDOCS_OUT_DIR}/reference_docs)
  file(MAKE_DIRECTORY ${READTHEDOCS_OUT_DIR}/reference_docs/reference_images)

  #install the python virtual environment with sphinx
  find_package(Sphinx REQUIRED)

  if(NOT DEFINED FOUR_C_SPHINX_THEME)
    set(FOUR_C_SPHINX_THEME sphinx_rtd_theme)
  endif()

  # configure readthedocs configuration
  configure_file(${READTHEDOCS_CONFIG_IN} ${READTHEDOCS_CONFIG_OUT} @ONLY)
  message("Sphinx configure file conf.in is written based on ${READTHEDOCS_CONFIG_IN}")

  add_custom_target(
    readthedocs ALL
    COMMAND create_rtd
    COMMAND ${baciname} -h > ${READTHEDOCS_OUT_DIR}/reference_docs/4C-help.txt
    COMMAND
      cmake ${PROJECT_SOURCE_DIR} --list-presets >
      ${READTHEDOCS_OUT_DIR}/reference_docs/4C-cmake-presets.txt
    COMMAND
      ${SPHINX_EXECUTABLE} -q -b html # -W turns warnings into errors!
      -i "${CMAKE_CURRENT_SOURCE_DIR}" -i "${READTHEDOCS_OUT_DIR}/reference_docs" -i
      "${PROJECT_SOURCE_DIR}/tests/framework-test" -i "${READTHEDOCS_TPLDEPEND_DIR}" -o
      "${READTHEDOCS_HTML_DIR}" -c "${READTHEDOCS_OUT_DIR}" -s "${PROJECT_BINARY_DIR}/doc/tmp"
    DEPENDS create_rtd ${baciname}
    COMMENT "Build readthedocs documentation with Sphinx"
    )

  list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
endif()
