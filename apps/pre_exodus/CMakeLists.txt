# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(OBJS_PRE_EXODUS
    # cmake-format: sortable
    4C_pre_exodus.cpp
    4C_pre_exodus_readbc.cpp
    4C_pre_exodus_reader.cpp
    4C_pre_exodus_validate.cpp
    4C_pre_exodus_writedat.cpp
    )

add_executable(pre_exodus EXCLUDE_FROM_ALL ${OBJS_PRE_EXODUS})
set_target_properties(pre_exodus PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} post)
four_c_set_up_executable(pre_exodus)
