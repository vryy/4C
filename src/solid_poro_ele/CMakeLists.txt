baci_auto_define_module()

set(_solid_poro_ele_dependencies
    # cmake-format: sortable
    lib
    solid_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_solid_poro_ele_dependencies})
