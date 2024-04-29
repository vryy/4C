four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    global_legacy_module
    inpar
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
