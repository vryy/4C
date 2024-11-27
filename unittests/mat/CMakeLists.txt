# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_mat)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_druckerprager_test.cpp
    4C_elasthyper_service_test.cpp
    4C_electrode_test.cpp
    4C_inelastic_defgrad_factors_test.cpp
    4C_material_service_test.cpp
    4C_multiplicative_split_defgrad_elasthyper_service_test.cpp
    4C_multiplicative_split_defgrad_elasthyper_test.cpp
    4C_muscle_utils_test.cpp
    4C_stvenantkirchhoff_test.cpp
    4C_viscoplastic_laws_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
