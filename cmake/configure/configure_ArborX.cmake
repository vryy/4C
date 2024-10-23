# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

message(STATUS "Fetch content for ArborX")
# Unconditionally turn on MPI support inside ArborX
set(ARBORX_ENABLE_MPI "ON")
fetchcontent_declare(
  arborx
  GIT_REPOSITORY https://github.com/arborx/ArborX.git
  GIT_TAG 10440a4a5e42c90a847800c5664934f7042febc9 #v1.4
  )
fetchcontent_makeavailable(arborx)

four_c_add_external_dependency(four_c_all_enabled_external_dependencies ArborX::ArborX)
