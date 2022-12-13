set(OBJS_DRT_PARTICLE_WALL
    # cmake-format: sortable
    particle_wall.cpp
    particle_wall_datastate.cpp
    particle_wall_discretization_runtime_vtu_writer.cpp
    particle_wall_result_test.cpp
    )

baci_add_library(alg_particle_wall ${OBJS_DRT_PARTICLE_WALL})

set(LIBS_ALG_PARTICLE_WALL
    # cmake-format: sortable
    drt_inpar
    drt_lib
    drt_particleengine
    )

target_link_libraries(alg_particle_wall PUBLIC ${LIBS_ALG_PARTICLE_WALL})
