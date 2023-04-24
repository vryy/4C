set(OBJS_PARTICLE_WALL
    # cmake-format: sortable
    particle_wall.cpp
    particle_wall_datastate.cpp
    particle_wall_discretization_runtime_vtu_writer.cpp
    particle_wall_result_test.cpp
    )

baci_add_library(alg_particle_wall ${OBJS_PARTICLE_WALL})

set(LIBS_ALG_PARTICLE_WALL
    # cmake-format: sortable
    inpar
    lib
    particleengine
    )

baci_add_dependency(alg_particle_wall ${LIBS_ALG_PARTICLE_WALL})
