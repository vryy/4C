baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    binstrategy
    config
    core
    cut
    fiber
    global_data
    inpar
    io
    mat
    nurbs_discret
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
