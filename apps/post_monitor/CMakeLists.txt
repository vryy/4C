# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(OBJS_POST_MONITOR
    # cmake-format: sortable
    4C_post_monitor.cpp
    )

add_executable(post_monitor EXCLUDE_FROM_ALL ${OBJS_POST_MONITOR})
set_target_properties(post_monitor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
four_c_set_up_executable(post_monitor)
