set(OBJS_DRT_PARTICLE_ALGORITHM
    # cmake-format: sortable
    particle_algorithm.cpp
    particle_algorithm_utils.cpp
    particle_dirichlet_bc.cpp
    particle_gravity.cpp
    particle_initial_field.cpp
    particle_input_generator.cpp
    particle_result_test.cpp
    particle_sim.cpp
    particle_temperature_bc.cpp
    particle_timint.cpp
    particle_viscous_damping.cpp
    )

baci_add_library(alg_particle_algorithm ${OBJS_DRT_PARTICLE_ALGORITHM})

set(LIBS_ALG_PARTICLE_ALGORITHM
    # cmake-format: sortable
    alg_particle_interaction
    alg_particle_rigidbody
    alg_particle_wall
    drt_adapter
    drt_lib
    drt_particleengine
    )

target_link_libraries(alg_particle_algorithm PUBLIC ${LIBS_ALG_PARTICLE_ALGORITHM})
