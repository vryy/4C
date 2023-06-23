baci_auto_define_module()

set(_lubrication_ele_dependencies
    # cmake-format: sortable
    core
    geometry
    inpar
    lib
    linalg
    mat
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_lubrication_ele_dependencies})
