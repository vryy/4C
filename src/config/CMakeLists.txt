# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# This module is special since files need to be generated.

# Import the global property with external dependency configuration as a variable.
get_property(
  FOUR_C_DEFINES_FOR_EXTERNAL_DEPENDENCIES GLOBAL PROPERTY FOUR_C_DEFINES_FOR_EXTERNAL_DEPENDENCIES
  )

# check if the compiler supports the visibility attribute
check_cxx_compiler_flag("-fvisibility=hidden" HAVE_VISIBILITY_ATTRIBUTE)

# configure a header file to pass some of the CMake settings to the source code
configure_file("4C_config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/config/4C_config.hpp")

four_c_get_git_revision_information()
# configure revision file to pass git revision information to the source code
configure_file(
  "4C_config_revision.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/config/4C_config_revision.hpp"
  )

set(_headers
    ${CMAKE_CURRENT_BINARY_DIR}/config/4C_config.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/config/4C_config_revision.hpp
    )

# Use the automatic definition which will create the usual targets.
four_c_auto_define_module()

# make sure the config subdirectory is visible to all compilation units
target_include_directories(
  ${AUTO_DEFINED_MODULE_NAME}_deps INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/config>
  )

# Manually add the configured headers to the config target.
target_sources(
  ${AUTO_DEFINED_MODULE_NAME}_deps
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_BINARY_DIR}
            FILES
            ${_headers}
  )
