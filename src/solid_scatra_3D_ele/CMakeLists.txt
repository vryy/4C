baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    solid_3D_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
