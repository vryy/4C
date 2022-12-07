set(OBJS_DRT_PARTICLE_ALGORITHM
    particle_algorithm_utils.cpp
    particle_algorithm.cpp
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
target_link_libraries(alg_particle_algorithm PUBLIC drt_adapter drt_particleengine drt_lib)
