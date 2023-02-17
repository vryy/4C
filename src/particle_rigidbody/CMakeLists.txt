set(OBJS_PARTICLE_RIGIDBODY
    # cmake-format: sortable
    particle_rigidbody.cpp
    particle_rigidbody_affiliation_pairs.cpp
    particle_rigidbody_datastate.cpp
    particle_rigidbody_result_test.cpp
    particle_rigidbody_runtime_vtp_writer.cpp
    )

baci_add_library(alg_particle_rigidbody ${OBJS_PARTICLE_RIGIDBODY})

set(LIBS_ALG_PARTICLE_RIGIDBODY
    # cmake-format: sortable
    io
    lib
    particleengine
    )

baci_add_dependency(alg_particle_rigidbody ${LIBS_ALG_PARTICLE_RIGIDBODY})
