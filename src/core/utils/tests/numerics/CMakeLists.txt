# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_numerics)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_cubic_spline_interpolation_test.cpp
    4C_utils_local_integration_test.cpp
    4C_utils_local_newton_test.cpp
    4C_utils_local_numeric_methods_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
