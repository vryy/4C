# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(OBJS_CREATERTD
    4C_create_rtdfiles_base.cpp 4C_create_rtdfiles_wrapper.cpp 4C_create_rtdfiles_utils.cpp
    )
add_executable(create_rtd EXCLUDE_FROM_ALL ${OBJS_CREATERTD})
set_target_properties(create_rtd PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
four_c_set_up_executable(create_rtd)
add_custom_target(rtd DEPENDS create_rtd)
