baci_auto_define_module()

set(_particle_rigidbody_dependencies
    # cmake-format: sortable
    inpar
    io
    lib
    particle_engine
    particle_interaction
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_particle_rigidbody_dependencies})
