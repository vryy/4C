# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

function(four_c_auto_define_python_bindings_submodule_tests)
  # only add tests if
  if(NOT FOUR_C_ENABLE_PYTHON_BINDINGS)
    return()
  endif()

  set(options "")
  set(oneValueArgs MODULE)
  set(multiValueArgs "")
  cmake_parse_arguments(
    _parsed
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
  if(DEFINED _parsed_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "There are unparsed arguments: ${_parsed_UNPARSED_ARGUMENTS}")
  endif()

  if(_parsed_MODULE)
    set(_module_under_test ${_parsed_MODULE})
  else()
    message(FATAL_ERROR "MODULE not set")
  endif()

  add_test(
    NAME ${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME}.${_module_under_test}
    COMMAND
      ${CMAKE_COMMAND} -E env bash -c
      "${FOUR_C_PYTHON_VENV_BUILD}/bin/python -m pytest ${CMAKE_CURRENT_SOURCE_DIR}"
    )

  # Tell pytest where to find the py4C module
  set_tests_properties(
    ${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME}.${_module_under_test}
    PROPERTIES ENVIRONMENT
               "PYTHONPATH=$<TARGET_FILE_DIR:${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME}>"
               LABELS
               minimal
    )

endfunction()
