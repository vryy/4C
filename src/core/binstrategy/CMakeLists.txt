four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    core
    inpar
    mortar
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
