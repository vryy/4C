baci_auto_define_module()

set(_particle_interaction_dependencies
    # cmake-format: sortable
    core
    inpar
    io
    lib
    linalg
    mat
    particle_algorithm
    particle_engine
    particle_wall
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_particle_interaction_dependencies})
