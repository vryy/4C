baci_auto_define_module()

set(_particle_algorithm_dependencies
    # cmake-format: sortable
    adapter
    inpar
    io
    lib
    particle_engine
    particle_interaction
    particle_rigidbody
    particle_wall
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_particle_algorithm_dependencies})
