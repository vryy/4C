baci_auto_define_module()

set(_solid_ele_dependencies
    # cmake-format: sortable
    lib
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_solid_ele_dependencies})
