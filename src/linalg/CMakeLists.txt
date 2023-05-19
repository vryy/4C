baci_auto_define_module()

set(_linalg_dependencies
    # cmake-format: sortable
    core
    lib
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_linalg_dependencies})
