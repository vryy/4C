baci_auto_define_module()

set(_discsh3_dependencies
    # cmake-format: sortable
    core
    inpar
    lib
    mat
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_discsh3_dependencies})
