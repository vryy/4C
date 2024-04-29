four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    fiber
    global_data
    inpar
    io
    lib
    mat
    nurbs_discret
    so3
    structure_new
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
