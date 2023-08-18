baci_auto_define_module()

set(_lubrication_ele_dependencies
    # cmake-format: sortable
    core
    inpar
    lib
    mat
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_lubrication_ele_dependencies})
