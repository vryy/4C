baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    core
    global_data
    inpar
    io
    lib
    mat
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
