four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    contact
    core
    global_data
    inpar
    io
    lib
    lubrication
    lubrication_ele
    mat
    mortar
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
