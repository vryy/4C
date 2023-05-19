baci_auto_define_module()

set(_particle_engine_dependencies
    # cmake-format: sortable
    binstrategy
    inpar
    io
    lib
    linalg
    particle_algorithm
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_particle_engine_dependencies})
