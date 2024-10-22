# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_particle_interaction)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_particle_interaction_dem_contact_normal_test.cpp
    4C_particle_interaction_dem_contact_rolling_test.cpp
    4C_particle_interaction_dem_contact_tangential_test.cpp
    4C_particle_interaction_sph_artificialviscosity_test.cpp
    4C_particle_interaction_sph_density_correction_test.cpp
    4C_particle_interaction_sph_equationofstate_test.cpp
    4C_particle_interaction_sph_kernel_test.cpp
    4C_particle_interaction_sph_momentum_formulation_test.cpp
    4C_particle_interaction_utils_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
