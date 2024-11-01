# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

message(STATUS "Fetch content for yaml-cpp")
fetchcontent_declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG f7320141120f720aecc4c32be25586e7da9eb978 # version 0.8.0
  )
fetchcontent_makeavailable(yaml-cpp)

four_c_add_external_dependency(four_c_all_enabled_external_dependencies yaml-cpp::yaml-cpp)
