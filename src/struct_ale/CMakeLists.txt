baci_auto_define_module()

set(_struct_ale_dependencies
    # cmake-format: sortable
    lib
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_struct_ale_dependencies})
