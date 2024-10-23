# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# This function links all default external dependencies to the given target.
# The keyword defines the usual linking visibility (PUBLIC, INTERFACE, PRIVATE)
function(four_c_link_default_external_libraries target keyword)
  target_link_libraries(${target} ${keyword} four_c_all_enabled_external_dependencies)
endfunction()
