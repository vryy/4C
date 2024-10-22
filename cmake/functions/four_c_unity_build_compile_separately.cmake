# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Exclude a file from unity build for its target.
function(four_c_unity_build_compile_separately _target _file)
  set_source_files_properties(
    ${_file} TARGET_DIRECTORY ${_target}_objs PROPERTIES SKIP_UNITY_BUILD_INCLUSION 1
    )
endfunction()
