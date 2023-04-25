set(OBJS_PARTICLE_ENGINE
    # cmake-format: sortable
    particle_engine.cpp
    particle_engine_communication_utils.cpp
    particle_engine_container.cpp
    particle_engine_container_bundle.cpp
    particle_engine_runtime_vtp_writer.cpp
    particle_engine_unique_global_id.cpp
    )

baci_add_library(alg_particle_engine ${OBJS_PARTICLE_ENGINE})

set(LIBS_ALG_PARTICLE_ENGINE
    # cmake-format: sortable
    alg_particle_algorithm
    io
    lib
    linalg
    )

baci_add_dependency(alg_particle_engine ${LIBS_ALG_PARTICLE_ENGINE})

set(OBJS_PARTICLE_ENGINE_LIB
    # cmake-format: sortable
    particle_engine_enums.cpp
    particle_engine_object.cpp
    )

baci_add_library(particleengine ${OBJS_PARTICLE_ENGINE_LIB})

set(LIBS_PARTICLEENGINE
    # cmake-format: sortable
    lib
    linalg
    )

baci_add_dependency(particleengine ${LIBS_PARTICLEENGINE})
