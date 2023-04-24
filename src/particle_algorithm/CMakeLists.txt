set(OBJS_PARTICLE_ALGORITHM
    # cmake-format: sortable
    particle_algorithm.cpp
    particle_algorithm_dirichlet_bc.cpp
    particle_algorithm_gravity.cpp
    particle_algorithm_initial_field.cpp
    particle_algorithm_input_generator.cpp
    particle_algorithm_result_test.cpp
    particle_algorithm_sim.cpp
    particle_algorithm_temperature_bc.cpp
    particle_algorithm_timint.cpp
    particle_algorithm_utils.cpp
    particle_algorithm_viscous_damping.cpp
    )

baci_add_library(alg_particle_algorithm ${OBJS_PARTICLE_ALGORITHM})

set(LIBS_ALG_PARTICLE_ALGORITHM
    # cmake-format: sortable
    adapter
    alg_particle_interaction
    alg_particle_rigidbody
    alg_particle_wall
    lib
    particleengine
    )

baci_add_dependency(alg_particle_algorithm ${LIBS_ALG_PARTICLE_ALGORITHM})
