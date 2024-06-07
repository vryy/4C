four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    binstrategy
    config
    core
    cut
    global_data
    inpar
    mat
    nurbs_discret
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
