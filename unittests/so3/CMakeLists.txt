# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_so3)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_so_contact_element_reference_configuration_test.cpp
    4C_so_element_service_test.cpp
    4C_so_hex8_determinant_analysis_test.cpp
    4C_so_hex8_test.cpp
    4C_so_tet4_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
