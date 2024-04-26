four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    binstrategy
    config
    core
    global_data
    inpar
    io
    lib
    particle_engine
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
