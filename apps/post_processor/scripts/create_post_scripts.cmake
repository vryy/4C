# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

function(copy_script script_name)
  add_custom_command(
    TARGET post_processor
    POST_BUILD
    COMMAND
      ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${script_name}
      ${PROJECT_BINARY_DIR}
    COMMENT "Create script ${script_name}"
    )
endfunction()

copy_script(post_ensight)
copy_script(post_gid)
copy_script(post_vti)
copy_script(post_vtu)
