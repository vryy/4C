baci_auto_define_module()

set(_lubrication_dependencies
    # cmake-format: sortable
    adapter
    inpar
    io
    lib
    lubrication_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_lubrication_dependencies})
