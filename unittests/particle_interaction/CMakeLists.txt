set(TESTNAME unittests_particle_interaction)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_particle_interaction_dem_contact_normal.cpp
    unit_particle_interaction_dem_contact_rolling.cpp
    unit_particle_interaction_dem_contact_tangential.cpp
    unit_particle_interaction_sph_artificialviscosity.cpp
    unit_particle_interaction_sph_density_correction.cpp
    unit_particle_interaction_sph_equationofstate.cpp
    unit_particle_interaction_sph_kernel.cpp
    unit_particle_interaction_sph_momentum_formulation.cpp
    unit_particle_interaction_utils.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} alg_particle_interaction)
