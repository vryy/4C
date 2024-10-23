# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(OBJS_PRE_LOCSYS
    # cmake-format: sortable
    4C_pre_locsys.cpp
    )

add_executable(pre_locsys EXCLUDE_FROM_ALL ${OBJS_PRE_LOCSYS})
set_target_properties(pre_locsys PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
four_c_set_up_executable(pre_locsys)
