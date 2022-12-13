set(OBJS_DRT_PARTICLE_ENGINE
    # cmake-format: sortable
    particle_communication_utils.cpp
    particle_container.cpp
    particle_container_bundle.cpp
    particle_engine.cpp
    particle_runtime_vtp_writer.cpp
    particle_unique_global_id.cpp
    )

baci_add_library(alg_particle_engine ${OBJS_DRT_PARTICLE_ENGINE})

set(LIBS_ALG_PARTICLE_ENGINE
    # cmake-format: sortable
    alg_particle_algorithm
    drt_io
    drt_lib
    linalg
    )

target_link_libraries(alg_particle_engine PUBLIC ${LIBS_ALG_PARTICLE_ENGINE})

set(OBJS_DRT_PARTICLE_ENGINE_LIB
    # cmake-format: sortable
    particle_enums.cpp
    particle_object.cpp
    )

baci_add_library(drt_particleengine ${OBJS_DRT_PARTICLE_ENGINE_LIB})

set(LIBS_DRT_PARTICLEENGINE
    # cmake-format: sortable
    drt_lib
    linalg
    )

target_link_libraries(drt_particleengine PUBLIC ${LIBS_DRT_PARTICLEENGINE})
