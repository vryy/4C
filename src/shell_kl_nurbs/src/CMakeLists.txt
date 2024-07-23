four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    global_data
    inpar
    mat
    structure_new
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
