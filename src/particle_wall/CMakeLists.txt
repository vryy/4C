baci_auto_define_module()

set(_particle_wall_dependencies
    # cmake-format: sortable
    binstrategy
    inpar
    io
    lib
    linalg
    particle_engine
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_particle_wall_dependencies})
