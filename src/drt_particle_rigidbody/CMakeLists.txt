set(OBJS_DRT_PARTICLE_RIGIDBODY
    particle_rigidbody_affiliation_pairs.cpp
    particle_rigidbody_datastate.cpp
    particle_rigidbody_result_test.cpp
    particle_rigidbody_runtime_vtp_writer.cpp
    particle_rigidbody.cpp
    )

baci_add_library(alg_particle_rigidbody ${OBJS_DRT_PARTICLE_RIGIDBODY})
target_link_libraries(alg_particle_rigidbody PUBLIC drt_lib drt_io drt_particleengine)
