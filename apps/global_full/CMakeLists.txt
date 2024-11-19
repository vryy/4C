# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(OBJS_FOUR_C_MAIN
    # cmake-format: sortable
    4C_global_full_cal_control.cpp
    4C_global_full_control.cpp
    4C_global_full_init_control.cpp
    4C_global_full_inp_control.cpp
    4C_global_full_main.cpp
    )

add_executable(${FOUR_C_EXECUTABLE_NAME} ${OBJS_FOUR_C_MAIN})
set_target_properties(
  ${FOUR_C_EXECUTABLE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
four_c_set_up_executable(${FOUR_C_EXECUTABLE_NAME})

# Create a symbolic link for the old project name
string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(FOUR_C_EXECUTABLE_NAME_LEGACY baci-${BUILD_TYPE})

add_custom_command(
  TARGET ${FOUR_C_EXECUTABLE_NAME}
  POST_BUILD
  COMMAND
    ${CMAKE_COMMAND} -E create_symlink ${FOUR_C_EXECUTABLE_NAME}
    ${CMAKE_BINARY_DIR}/${FOUR_C_EXECUTABLE_NAME_LEGACY}
  COMMENT "Creating symbolic link"
  )

# Create a json schema file for the input
add_custom_command(
  TARGET ${FOUR_C_EXECUTABLE_NAME}
  POST_BUILD
  COMMAND
    ${FOUR_C_ENABLE_ADDRESS_SANITIZER_TEST_OPTIONS} "$<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}>" -p >
    ${CMAKE_BINARY_DIR}/4C_metadata.yaml
  COMMAND
    ${FOUR_C_PYTHON_VENV_BUILD}/bin/python ${PROJECT_SOURCE_DIR}/utilities/create_json_schema.py
    ${CMAKE_BINARY_DIR}/4C_metadata.yaml ${CMAKE_BINARY_DIR}/4C_schema.json
  COMMENT "Creating schema file for input"
  )
