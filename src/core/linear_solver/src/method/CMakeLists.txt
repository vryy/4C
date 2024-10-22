# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

four_c_auto_define_module()

# Clashing definitions of `View` inside Trilinos
four_c_unity_build_compile_separately(
  ${AUTO_DEFINED_MODULE_NAME} 4C_linear_solver_method_direct.cpp
  )
