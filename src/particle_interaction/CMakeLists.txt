baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    global_data
    inpar
    io
    lib
    mat
    particle_algorithm
    particle_engine
    particle_wall
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
