set(OBJS_DRT_PARTICLE_RIGIDBODY
    # cmake-format: sortable
    particle_rigidbody.cpp
    particle_rigidbody_affiliation_pairs.cpp
    particle_rigidbody_datastate.cpp
    particle_rigidbody_result_test.cpp
    particle_rigidbody_runtime_vtp_writer.cpp
    )

baci_add_library(alg_particle_rigidbody ${OBJS_DRT_PARTICLE_RIGIDBODY})

set(LIBS_ALG_PARTICLE_RIGIDBODY
    # cmake-format: sortable
    drt_io
    drt_lib
    drt_particleengine
    )

target_link_libraries(alg_particle_rigidbody PUBLIC ${LIBS_ALG_PARTICLE_RIGIDBODY})
