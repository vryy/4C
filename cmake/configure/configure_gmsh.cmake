# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# find gmsh (lower case for the package name)
find_package(gmsh REQUIRED)
if(TARGET gmsh::lib AND WIN32)
  set_property(
    TARGET gmsh::lib
    APPEND
    PROPERTY INTERFACE_LINK_LIBRARIES winmm
    )
  add_library(gmsh_lib ALIAS gmsh::lib)
elseif(TARGET gmsh::shared)
  add_library(gmsh_lib ALIAS gmsh::shared)
endif()

if(gmsh_FOUND)
  get_target_property(GMSH_INCLUDES gmsh_lib INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "gmsh include directories: ${GMSH_INCLUDES}")

  target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE gmsh_lib)

  four_c_remember_variable_for_install(FOUR_C_GMSH_ROOT)
endif()
