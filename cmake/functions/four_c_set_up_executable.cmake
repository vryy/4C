# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Call this function on an executable of this project.
# The executable will link to the main library and use the internal compiler settings.
function(four_c_set_up_executable target)
  target_link_libraries(${target} PRIVATE ${FOUR_C_LIBRARY_NAME})
  target_link_libraries(${target} PRIVATE four_c_private_compile_interface)
endfunction()
