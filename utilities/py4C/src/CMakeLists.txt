# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# create python bindings module
# NOTE: The NO_EXTRAS option is required for GCC to prevent it from adding the -flto flag, which would break pybind11
pybind11_add_module(${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME} NO_EXTRAS 4C_python_bindings_main.cpp)

# Add python bindings to the full target list
add_dependencies(full ${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME})

add_subdirectory(config)

set_target_properties(
  ${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME}
  PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY
    "${PROJECT_BINARY_DIR}/${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME}/${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME}"
    # Keep the binary tiny
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN ON
  )
install(TARGETS ${FOUR_C_PYTHON_BINDINGS_PROJECT_NAME} DESTINATION .)
