set(OBJS_DRT_PARTICLE_ENGINE
    particle_communication_utils.cpp
    particle_container.cpp
    particle_container_bundle.cpp
    particle_engine.cpp
    particle_runtime_vtp_writer.cpp
    particle_unique_global_id.cpp
    )

baci_add_library(alg_particle_engine ${OBJS_DRT_PARTICLE_ENGINE})
target_link_libraries(alg_particle_engine PUBLIC drt_lib linalg drt_io)

set(OBJS_DRT_PARTICLE_ENGINE_LIB particle_enums.cpp particle_object.cpp)
baci_add_library(drt_particleengine ${OBJS_DRT_PARTICLE_ENGINE_LIB})
target_link_libraries(drt_particleengine PUBLIC drt_lib)
